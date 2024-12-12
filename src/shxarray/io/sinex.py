import numpy as np
import xarray as xr
import shxarray
import os
import gzip
from shxarray.core.sh_indexing import SHindexBase
from datetime import datetime,timedelta
from shxarray.core.logging import logger

def sinex2date(snxdate:str):
    yr,doy,sec=[int(x) for x in snxdate.split(":")]
    return datetime(yr+2000 if yr<50 else yr+1900,1,1)+timedelta(days=doy-1,seconds=sec)



def read_vec(fileobj,dsout,blockname):
    svtype=None
    vtype=None
    if blockname == "SOLUTION/ESTIMATE":
        vtype='sol_est'
        svtype='sol_std'
    elif blockname == "SOLUTION/APRIORI":
        vtype='apri_est'

    elif blockname == 'SOLUTION/NORMAL_EQUATION_VECTOR':
        vtype='rhs'

    nest=dsout.dims['nm']
    
    ids=np.empty([nest],dtype=int)
    #allocate space
    nm=np.empty([nest],dtype=object)
    est=np.zeros([nest])
    if svtype:
        std=np.zeros([nest])
    epoch=np.empty([nest],dtype=datetime)
    i=0
    for line in fileobj:
        if line.startswith('-'):
            #end of block encountered
            break
        elif line.startswith('*'):
            #comment
            continue
        
        fields = line.split()
        ids[i]=int(fields[0])
        n = int(fields[2])
        est[i]=float(fields[8])
        if svtype:
            std[i]=float(fields[9])

        m=int(fields[4])
        if fields[1] == "SN":
            m = -m
        elif fields[1] != "CN":
            raise NotImplementedError(f"SINEX Parameter {fields[1]} is currently not supported")
        nm[i]=(n,m)
        epoch[i]=sinex2date(fields[5])
        i+=1
    #update dsout
    dsout[vtype]=('nm',est)
    if svtype:
        dsout[svtype]=('nm',std)
    
    # ignore parameter specific time epochs for now 
    # dsout[f'{vtype}_epoch']=('nm',epoch)


    if 'nm' not in dsout.indexes:
        #possibly build index (if it has not been build already)
        mi=SHindexBase.mi_fromtuples(nm)
        mi=xr.Coordinates.from_pandas_multiindex(mi, "nm")
        
        dsout=dsout.assign_coords(mi)
    
    return dsout

def read_symmat(file_or_obj,dsout,blockname):
    nest=dsout.dims['nm']
    mat=np.zeros([nest,nest],order='C')
    for line in file_or_obj:
        if line.startswith('-'):
            #end of block encountered
            break
        elif line.startswith('*'):
            #comment
            continue

        data=[float(x) for x in line.split()]

        irow=int(data[0])-1 #note zero indexing
        icol=int(data[1]) -1
        ndat=len(data)-2
        mat[irow,icol:icol+ndat]=data[2:]
    #mirror the upper triangle in the lower part
    mat=np.triu(mat,k=1).T+mat

    if "nm_" not in dsout.indexes:
        #add the transposed index
        mi_=SHindexBase.mi_toggle(dsout.indexes['nm'])
        mi_=xr.Coordinates.from_pandas_multiindex(mi_, "nm_")
        dsout=dsout.assign_coords(mi_)

    dsout['N']=(['nm','nm_'],mat)
    breakpoint()
    return dsout
# dictionary to lookup functions to dispatch the block parsing to (note some functions are the same for different blocks)
blockdispatch={"SOLUTION/ESTIMATE":read_vec,'SOLUTION/APRIORI':read_vec,
               'SOLUTION/NORMAL_EQUATION_VECTOR':read_vec,
               'SOLUTION/NORMAL_EQUATION_MATRIX U':read_symmat,
               'SOLUTION/NORMAL_EQUATION_MATRIX L':read_symmat}

compatversions=["2.02"]

def read_sinex(file_or_obj,stopatmat=False):
    needsClosing=False
    if type(file_or_obj) == str:
        needsClosing=True
        if file_or_obj.endswith('.gz'):
            file_or_obj=gzip.open(file_or_obj,'rt')
        else:
            file_or_obj=open(file_or_obj,'rt')
    
    # read first line
    header=file_or_obj.readline().split()
    if header[1] not in compatversions:
        raise RuntimeError(f"read_sinex is not compatible with {headerline[1]}")
    
    nest=int(header[-3])
    tstart=sinex2date(header[5])
    tend=sinex2date(header[6])
    #initialize xarray dataset with some scalar vars to augment
    dsout=xr.Dataset(dict(tstart=tstart,tend=tend,snx_ids=("nm",np.arange(1,nest+1))))
    
    

    # loop until a block is encountered and then dispatch to appropriate function
    for line in file_or_obj:
        if line[0] == "+":
            block=line[1:].strip()

            if block not in blockdispatch.keys():
                logger.info(f"Ignoring block {block}")
                continue
            if stopatmat and "MATRIX" in block:
                logger.info(f"Encountered {block}, stopping")
                break
            logger.info(f"Reading block {block}")
            dsout=blockdispatch[block](file_or_obj,dsout,block)
        
    if needsClosing:
        file_or_obj.close()
    
    return dsout 
    
