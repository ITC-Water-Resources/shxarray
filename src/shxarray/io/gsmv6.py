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
from itertools import chain
from datetime import datetime,timedelta


def read_yaml_header(buf):

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
    return attr

def read_cnes_header(buf):

    #setup global attributes
    attr={}
    buf.seek(0)
    for ln in buf:
        if re.search(b'FIRST.*CNES.*$',ln):
            #extract start and end time from the header line
            tst,tnd=ln.decode('utf-8').split("_")[1].split("-")
            attr['tstart']=datetime(int(tst[0:4]),1,1)+timedelta(days=int(tst[4:8])-1)
            attr['tend']=datetime(int(tnd[0:4]),1,1)+timedelta(days=int(tnd[4:8])-1)
            attr['source']=ln[6:].decode('utf-8')
        elif re.search(b'^SHM',ln):
            attr['nmax']=int(ln.split()[1])
        elif re.search(b'^EARTH',ln):
            gm,rearth=[float(val) for val in ln.split()[1:]]
            attr["gm"]=gm
            attr["re"]=rearth
            
    return attr

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
    headertype=None
    for ln in fileobj:
        if b'# End of YAML header' in ln:
            lineiterator=fileobj
            headertype='YAML'
            break
        elif b'GRCOF2' in ln:
            #old style GSM header
            #include first valid line in the line iterator
            lineiterator=chain([ln],fileobj)
            headertype='CNES'
            break
        else:
            buf.write(ln)

    if headertype == "YAML":
        attr=read_yaml_header(buf)
    elif headertype=='CNES':
        attr=read_cnes_header(buf)
    else:
        raise IOError("header not recognized")

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
    for ln in lineiterator:
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
    
    nm_mi=SHindexBase.mi_fromtuples(nm)
    if hasattr(xr,'Coordinates'):
        #only in newer xarray versions: convert to xarray coordinates
        coords=xr.Coordinates.from_pandas_multiindex(nm_mi, SHindexBase.name)
        if time:
            coords.update(dict(time=time))
    else:
        coords={SHindexBase.name:nm_mi}
        if time:
            coords["time"]=time

    if time:
        shp=["time",SHindexBase.name]

        # coords={SHindexBase.name:nm_mi,"time":time}
        #also expand variables
        cnm=np.expand_dims(cnm[0:ncount], axis=0)
        sigcnm=np.expand_dims(sigcnm[0:ncount],axis=0)
    else:
        shp=[SHindexBase.name]
        # coords={SHindexBase.name:nm_mi}
        cnm=cnm[0:ncount]
        sigcnm=sigcnm[0:ncount]
    
    ds=xr.Dataset(data_vars=dict(cnm=(shp,cnm,get_cfatts("stokes")),sigcnm=(shp,sigcnm,get_cfatts("stokes stdv"))),coords=coords,attrs=attr)
    return ds
