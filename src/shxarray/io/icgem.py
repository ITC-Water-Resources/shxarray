# This file is part of frommle2.
# frommle2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# frommle2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2021

import gzip
import xarray as xr
import re
import sys
import numpy as np
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.logging import logger 
from datetime import datetime,timedelta
from shxarray.core.cf import get_cfatts

def readIcgem(fileobj,nmaxstop=sys.maxsize):
    needsClosing=False
    if type(fileobj) == str:
        needsClosing=True
        if fileobj.endswith('.gz'):
            fileobj=gzip.open(fileobj,'rb')
        else:
            fileobj=open(fileobj,'rb')

    #first read the icgem header 
    inheader=False
    hdr={}
    for ln in fileobj:
        if b'begin_of_head' in ln:
            inheader=True
            continue
        if b'end_of_head' in ln:
            break
        
        spl=ln.decode('utf-8').split()
        if len(spl) == 2:
            #insert name value pairs in the hdr dict
            hdr[spl[0]]=spl[1]
    
    #extract relevant parameters from the header
    attr={}
    try:
        nmaxsupp=int(hdr["max_degree"])
        attr["nmaxfile"]=nmaxsupp
        if nmaxsupp < nmaxstop:
            attr["nmax"]=nmaxsupp
        else:
            attr["nmax"]=nmaxstop
        nmax=attr["nmax"]

        if nmax > nmaxsupp:
            logger.warning("Nmax ({nmax}) requested larger than supported, higher degree coefficients will be set to zero")


        if 'format' in hdr:
            attr["format"]=hdr['format']
        else:
            attr["format"]="icgem"
        
        if "norm" in hdr:
            attr["norm"]=hdr["norm"]
        
        attr["gm"]=float(hdr["earth_gravity_constant"])
        attr["re"]=float(hdr["radius"])
        attr["modelname"]=hdr["modelname"]
    except KeyError:
    #some values may not be present but that is ok
        pass
    
    #Non standard HACK to try to retrieve the epoch from the modelname (GRAZ monthly solutions only)
    if "modelname" in hdr:
        try:
            time=[datetime.strptime(hdr['modelname'][-7:],"%Y-%m")+timedelta(days=14.5)]
        except ValueError:
            time=None
    else:
        time=None


    nsh=SHindexBase.nsh(nmax,squeeze=True)
    cnm=np.zeros([nsh])
    sigcnm=np.zeros([nsh])
    ncount=0
    nm=[]
    #continue reading the data
    dataregex=re.compile(b'^gfc')
    for ln in fileobj:
        if dataregex.match(ln):
            lnspl=ln.replace(b"D",b"E").split()
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
