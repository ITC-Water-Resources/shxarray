import numpy as np
import xarray as xr
import shxarray
import os
import gzip
from shxarray.core.sh_indexing import SHindexBase
from datetime import datetime,timedelta
from shxarray.core.logging import shxlogger
from shxarray.io.gzipwrap import gzip_open_r
def sinex2date(snxdate:str)->datetime:
    """
        Convert sinex datestring in yy:doy:seconds to python datetime

    Parameters
    ----------
    snxdate : str
        datestring in sinex format
        

    Returns
    -------
    datetime
       Specified date and time as datetime object 

    """
    yr,doy,sec=[int(x) for x in snxdate.split(":")]
    return datetime(yr+2000 if yr<50 else yr+1900,1,1)+timedelta(days=doy-1,seconds=sec)



def read_vec(fileobj,dsout,blockname):
    """
        Read a SINEX block containing a vector

    Parameters
    ----------
    fileobj : 
        io buffer to read lines from
        
    dsout : xarray.Dataset
        xarray.Dataset to augment the vector data to
        
        
    blockname : str
        name of the SINEX block. should be one of:
        SOLUTION/ESTIMATE
        SOLUTION/APRIORI
        SOLUTION/NORMAL_EQUATION_VECTOR
        

    Returns
    -------
    an updated xarray.Dataset holding the new data in new variables:
    depending on the blockname these ar: sol_est,sol_std,apri_est or rhs

    """
    svtype=None
    vtype=None
    if blockname == "SOLUTION/ESTIMATE":
        vtype='sol_est'
        svtype='sol_std'
    elif blockname == "SOLUTION/APRIORI":
        vtype='apri_est'

    elif blockname == 'SOLUTION/NORMAL_EQUATION_VECTOR':
        vtype='rhs'
    else:
        raise NotImplementedError(f"Cannot process {blockname}")
    nest=dsout.sizes['nm']
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
        dsout=dsout.sh.set_nmindex(mi)
    
    return dsout

def read_symmat(fileobj,dsout,blockname):
    """
        Reads a triangular matrix from a SINEX block and returns a symmetric version
    Parameters
    ----------
    fileobj : 
        io buffer to read lines from
        
    dsout : xarray.Dataset
        xarray.Dataset to augment the matrix data to
        
        
    blockname : str
        name of the SINEX block. should be one of:
        SOLUTION/NORMAL_EQUATION_MATRIX U
        SOLUTION/NORMAL_EQUATION_MATRIX L
        

    Returns
    -------
    an updated xarray.Dataset holding the new matrix in a new variable 'N'

    """
    if not blockname.startswith('SOLUTION/NORMAL_EQUATION_MATRIX'):
        raise RuntimeError(f"Wrong block {blockname}?")
    nest=dsout.sizes['nm']
    mat=np.zeros([nest,nest],order='C')
    data=np.zeros([3])
    for line in fileobj:
        if line[0] == '-':
            #end of block encountered
            break
        elif line[0] =='*':
            #comment
            continue
        data=[float(x) for x in line.split()]
        irow=int(data[0])-1 #note zero indexing
        icol=int(data[1])-1
        ndat=len(data)-2
        mat[irow,icol:icol+ndat]=data[2:]
    #mirror the upper triangle in the lower part
    mat=np.triu(mat,k=1).T+mat

    if "nm_" not in dsout.indexes:
        #add the transposed index
        mi_=SHindexBase.mi_toggle(dsout.indexes['nm'])
        dsout=dsout.sh.set_nmindex(mi_,'_')

    dsout['N']=(['nm','nm_'],mat)
    return dsout

def read_statistics(fileobj,dsout,blockname):
    """
        Reads a SINEX block holding statistical metrics
    Parameters
    ----------
    fileobj : 
        io buffer to read lines from
        
    dsout : xarray.Dataset
        xarray.Dataset to augment the matrix data to
        
        
    blockname : str
        name of the SINEX block. should be SOLUTION/STATISTICS
        

    Returns
    -------
    an updated xarray.Dataset holding the available statistics  as scalar variables


    """
    if blockname != 'SOLUTION/STATISTICS':
        raise RuntimeError(f"Wrong block encountered?")

    for line in fileobj:
        if line.startswith('-'):
            #end of block encountered
            break
        elif line.startswith('*'):
            #comment
            continue
            continue
        
        if line.startswith(" NUMBER OF DEGREES OF FREEDOM"):
            varname="dof"
            tp=int
        elif line.startswith(" NUMBER OF OBSERVATIONS"):
            varname="nobs"
            tp=int
        elif line.startswith(" NUMBER OF UNKNOWNS"):
            varname="nunknown"
            tp=int
        elif line.startswith(" WEIGHTED SQUARE SUM OF O-C"):
            varname="ltpl"
            tp=float
        elif line.startswith(" VARIANCE FACTOR"):
            varname="sigma0_2"
            tp=float
        elif line.startswith(" SQUARE SUM OF RESIDUALS (VTPV)"):
            varname="ltpl"
            tp=float
        else:
            tp=None
            varname =None
            shxlogger.warning(f"ignoring {blockname} entry {line}")

        if varname:
            spl = line.split()
            val=tp(spl[-1])
            dsout[varname]=val     
    return dsout




# dictionary to lookup functions to dispatch the block parsing to (note some functions are the same for different blocks)
blockdispatch={"SOLUTION/ESTIMATE":read_vec,'SOLUTION/APRIORI':read_vec,
               'SOLUTION/NORMAL_EQUATION_VECTOR':read_vec,
               'SOLUTION/STATISTICS':read_statistics,
               'SOLUTION/NORMAL_EQUATION_MATRIX U':read_symmat,
               'SOLUTION/NORMAL_EQUATION_MATRIX L':read_symmat}

compatversions=["2.02"]

def read_sinex(file_or_obj,stopatmat=False):
    """
        Reads normal equation information from a SINEX file
    Parameters
    ----------
    file_or_obj : 
        IO buffer or filename with the SINEX data source
        
    stopatmat : bool
        Stop reading from the source when encountering a MATRIX block to speed up. (default=False)
        

    Returns
    -------
       a xarray.Dataset holding the normal equation information 
        

    """
    needsClosing=False
    if type(file_or_obj) == str:
        needsClosing=True
        if file_or_obj.endswith('.gz'):
            file_or_obj=gzip_open_r(file_or_obj,textmode=True)
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
                shxlogger.info(f"Ignoring block {block}")
                continue
            if stopatmat and "MATRIX" in block:
                shxlogger.info(f"Encountered {block}, stopping")
                break
            shxlogger.info(f"Reading block {block}")
            dsout=blockdispatch[block](file_or_obj,dsout,block)
        
    if needsClosing:
        file_or_obj.close()
    
    return dsout 
    
