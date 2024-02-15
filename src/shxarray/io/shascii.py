# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

import xarray as xr
from shxarray.core.time import decyear2dt
import numpy as np
import sys
import gzip
from shxarray.core.sh_indexing import SHindexBase
import re

# def writeSHAscii(fileobj,ds,cnmv='cnm',sigcnmv=None):
    # """Writes a dataset array with sh data to an ascii file"""
    
    # needsClosing=False

    # if type(fileobj) == str:
        # needsClosing=True
        # if fileobj.endswith('.gz'):
            # fileobj=gzip.open(fileobj,'wt')
        # else:
            # fileobj=open(fileobj,'wt')
    # nmax=ds.sh.nmax

    # #TODO extract time epochs in decimal years from the data
    # tstart=0.0
    # tcent=0.0
    # tend=0.0

    # fileobj.write(f" META    {nmax}    {tstart}   {tcent}   {tend}\n") 
    # #loop over all coefficients (make sure to sort them appropriately)
    # sortds=ds.sortby(['n','m','t'])
    # cnmvals=np.zeros([2])
    # ncoef=0
    # for idx,el in zip(sortds.shg,sortds[cnmv].values):
        # n,m,t=idx.data[()]
        # cnmvals[t]=el

        # if m == 0:
            # #no need to wait for a sine coefficient
            # ncoef=2
            # cnmvals[1]=0.0
        # else:
            # ncoef+=1
        
        # if ncoef == 2: 
            # fileobj.write(f"{n:6d} {m:6d} {cnmvals[0]:17.10e} {cnmvals[1]:17.10e}\n")
            # ncoef=0
            # cnmvals[:]=0.0

    # if needsClosing:
        # fileobj.close()


def readSHAscii(fileobj,nmaxstop=sys.maxsize):

    if type(fileobj) == str:
        if fileobj.endswith(".gz"):
            fid=gzip.open(fileobj,'rt')
        else:
            fid=open(fileobj,'rt')
    else:
        fid=fileobj
    
    #defaults
    attr={}
    nmaxfile=-1 # unknown
    
    #Try seeing if a META header line is present 
    frstln=fid.readline()
    if re.search(' +META',frstln):
        metaln=frstln.split()
        nmaxfile=int(metaln[1])

        if nmaxstop <= nmaxfile:
            nmax=nmaxstop
        else:
            nmax=nmaxfile

        #Extract time tags
        stime,ctime,etime=[decyear2dt(float(x)) for x in metaln[2:5]]
        attr={"nmax":nmax,"nmaxfile":nmaxfile,"CTime":ctime,"STime":stime,"ETime":etime}
        #read next line already (may already contain data)
        frstln=fid.readline()
   
    #check how many columns are present 
    ln=frstln.split()
    nd=len(ln)
    if nd == 4:
        sigma=False
    elif nd == 6:
        sigma=True
    else:
        raise IOError(f"Unexpected amount of columns: {len(ln)}")
    if nmaxfile > 0: 
        nshmax=SHindexBase.nsh(nmaxfile)
    else:
        nshmax=sys.maxsize

    ncount=0 
    n,m=[int(x) for x in ln[0:2]]
    coefs=[float(x) for x in ln[2:nd]]
    
    # Add cosine entry
    nm=[(n,m)]
    cnm=[coefs[0]]
    if nd == 6:
        sigcnm=[coefs[2]]
    ncount+=1 
    # add Sine entry if order is not zero    
    if m > 0:
        nm.append((n,-m))
        cnm.append(coefs[1])
        if nd == 6:
            sigcnm.append(coefs[3])
        ncount+=1 
    

    for lnstr in fid:
        ln=lnstr.split()
        if len(ln) == 0:
            break
        n,m=[int(x) for x in ln[0:2]]
        coefs=[float(x) for x in ln[2:nd]]
    
        # Add cosine entry
        nm.append((n,m))
        cnm.append(coefs[0])
        if nd == 6:
            sigcnm.append(coefs[2])
        ncount+=1 
        # add Sine entry if order is not zero    
        if m > 0:
            nm.append((n,-m))
            cnm.append(coefs[1])
            if nd == 6:
                sigcnm.append(coefs[3])
            ncount+=1 
        if n > nmaxstop and ncount > nshmax:
            #all requested data has been read already
            break
         
    #create an xarray dataset
    
    if nd == 6:
        ds=xr.Dataset(data_vars=dict(cnm=([SHindexBase.name],cnm),sigcnm=([SHindexBase.name],sigcnm)),coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm)},attrs=attr)
    else:
        ds=xr.Dataset(data_vars=dict(cnm=([SHindexBase.name],cnm)),coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm)},attrs=attr)
    return ds
