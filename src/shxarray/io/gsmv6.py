# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import gzip
import xarray as xr
from shxarray.core.sh_indexing import SHindexBase
import re
import sys
from io import BytesIO
import yaml
import numpy as np
from shxarray.core.cf import get_cfatts

def readGSMv6(fileobj,nmaxstop=sys.maxsize):
    needsClosing=False
    if type(fileobj) == str:
        needsClosing=True
        if fileobj.endswith('.gz'):
            fileobj=gzip.open(fileobj,'rb')
        else:
            fileobj=open(fileobj,'rb')

    #first read the yaml header 
    buf=BytesIO()
    for ln in fileobj:
        if b'# End of YAML header' in ln:
            break
        else:
            buf.write(ln)

    hdr=yaml.safe_load(buf.getvalue())["header"]
    
    #setup global attributes
    attr={}
    attr["nmaxfile"]=hdr["dimensions"]["degree"]
    if not "nmax" in attr:
        attr["nmax"]=attr["nmaxfile"]

    attr["tstart"]=hdr["global_attributes"]["time_coverage_start"]
    attr["tend"]=hdr["global_attributes"]["time_coverage_end"]

    nonstand=hdr["non-standard_attributes"]

    attr["gm"]=nonstand["earth_gravity_param"]["value"]
    attr["re"]=nonstand["mean_equator_radius"]["value"]
    

    nmax=attr["nmax"]


    nsh=SHindexBase.nsh(nmax,squeeze=True)
    
    cnm=np.zeros([nsh])
    sigcnm=np.zeros([nsh])
    if "tstart" in attr and "tend" in attr:
        time=[attr['tstart']+(attr['tend']-attr['tstart'])/2]
    else:
        time=None

    ncount=0
    nm=[]
    #continue reading the data
    dataregex=re.compile(b'^GRCOF2')
    for ln in fileobj:
        if dataregex.match(ln):
            lnspl=ln.split()
            n=int(lnspl[1])
            if n> nmaxstop:
                if ncount > nsh:
                    #all required coefficients have been read (no need to read the file further)
                    break
                continue
                

            m=int(lnspl[2])

            cnm[ncount]=float(lnspl[3])
            sigcnm[ncount]=float(lnspl[5])
            nm.append((n,m))
            ncount+=1
            
            #possibly also add snm coefficients
            if m!=0:
                cnm[ncount]=float(lnspl[4])
                sigcnm[ncount]=float(lnspl[6])
                nm.append((n,-m))
                ncount+=1

    if needsClosing:
        fileobj.close()
    
    if time:
        shp=["time",SHindexBase.name]
        coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm),"time":time}
        #also expand variables
        cnm=np.expand_dims(cnm[0:ncount], axis=0)
        sigcnm=np.expand_dims(sigcnm[0:ncount],axis=0)
    else:
        shp=[SHindexBase.name]
        coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm)}
        cnm=cnm[0:ncount]
        sigcnm=sigcnm[0:ncount]
    
    ds=xr.Dataset(data_vars=dict(cnm=(shp,cnm,get_cfatts("stokes")),sigcnm=(shp,sigcnm,get_cfatts("stokes stdv"))),coords=coords,attrs=attr)
    return ds
