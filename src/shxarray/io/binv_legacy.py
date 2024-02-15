# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#






import struct
import numpy as np
import re
from shxarray.core.sh_indexing import SHindexBase
import xarray as xr
import sparse

def getBDcoords(ddict,trans):

    if not ddict['type'] in ["BDFULLV0","BDFULLVN"]:
        raise RuntimeError("cannot get coordinates for this type of matrix (yet)")

    mxblk=np.max(np.diff(ddict['blockind']))
    if trans:
        mshx,mshy=np.meshgrid([x for x in range(mxblk)],[y for y in range(mxblk)])
    else:
        mshy,mshx=np.meshgrid([x for x in range(mxblk)],[y for y in range(mxblk)])
    # mshx=mshx.reshape([mshx.size])
    # mshy=mshy.reshape([mshy.size])

    #initialize coords array
    coords=np.zeros([2,ddict['pval1']*ddict['pval2']],dtype=np.int64)

    blkstrt=0
    nshift=0
    for blkend in ddict['blockind']:
        blklen=(blkend-blkstrt)
        sz=blklen**2
        coords[0,nshift:nshift+sz]=mshx[0:blklen,0:blklen].reshape([sz])+blkstrt
        coords[1,nshift:nshift+sz]=mshy[0:blklen,0:blklen].reshape([sz])+blkstrt
        # import pdb;pdb.set_trace()
        #update blockstart
        blkstrt=blkend
        nshift+=sz
    return coords

def get_shmi(charar):
    shregex=re.compile('^G[CS]N')
    nm=[]
    for el in charar:
        if shregex.search(el):
            n=int(el[4:7])
            m=int(el[7:10])
            if el[1:2] == 'S':
                m=-m
        else:
            raise KeyError("Error parsing SH coefficient")
        nm.append((n,m))
    nmmi=SHindexBase.mi_fromtuples(nm)
    return nmmi

def readBINV(file_or_obj,trans=False,nmax=-1):
    """Reads in a legacy binary file written using the fortran RLFTlbx"""
    dictout={}
    #default to assuming the file is in little endian
    endianness='<'
    if type(file_or_obj) == str:
        #open filename in binary mode
        fid=open(file_or_obj,'rb')
    else:
        fid=file_or_obj
    #read the endianess checker
    (endian,)=struct.unpack(endianness+"H", fid.read(2))
    #compare the magic number (should be 18754 on a system with the same endiannes of the file)
    if endian != 18754:
        #switch to big endian
        endianness='>'

    dictout['version']='BI'+fid.read(6).decode('utf-8')
    vnum=float(dictout['version'][4:7])

    #read type, description etc

    dictout['type']=fid.read(8).decode('utf-8')
    dictout['description']=fid.read(80).decode('utf-8')

    #read integer meta information
    (nints,ndbls,nval1,nval2)=struct.unpack(endianness+'IIII',fid.read(4*4))

    if vnum < 2.4:
        (pval1,pval2)=struct.unpack(endianness+'II',fid.read(4*2))
    else:
        (pval1,pval2)=struct.unpack(endianness+'LL',fid.read(8*2))


    if vnum <= 2.1:
        if dictout["type"] in ['SYMV0___','BDFULLV0','BDSYMV0','BDFULLVN']:
            nvec=0
            pval2=1
        elif dictout["type"] == "SYMV1___":
            nvec=1
            pval2=1
        elif dictout["type"] == "SYMV2___":
            nvec=2
            pval2=1
        elif dictout["type"] == "FULLSQV0":
            nvec=0
            pval2=pval1

        nread=0
        nval2=nval1
    else:
        (nvec,nread)=struct.unpack(endianness+"II",fid.read(4*2))

    # dictout["nval1"]=nval1
    # dictout["nval2"]=nval2
    dictout["pval1"]=pval1
    dictout["pval2"]=pval2

    #read type dependent index data
    if dictout['type'] in ["BDSYMV0_","BDSYMVN_","BDFULLV0","BDFULLVN"]:
        (nblocks,)=struct.unpack(endianness+'I',fid.read(4))
        dictout["nblocks"]=nblocks

    if nread >0:
        dictout["readme"]=fid.read(nread*80).decode('utf-8')

    names=[]
    vals=[]
    if nints > 0:
        inames=np.fromfile(fid,dtype='|S24',count=nints).astype('|U24')
        names.extend(list(inames))
        if vnum <= 2.4:
            ivals=np.fromfile(fid,dtype=endianness+'I',count=nints)
        else:
            ivals=np.fromfile(fid,dtype=endianness+'L',count=nints)
        vals.extend(ivals)
    if ndbls >0:
        dnames=np.fromfile(fid,dtype='|S24',count=ndbls).astype('|U24')
        names.extend(list(dnames))
        dvals=np.fromfile(fid,dtype=endianness+'d',count=ndbls)
        vals.extend(dvals)

    if names:
        dictout.update({ky.strip():val for ky,val in zip(names,vals)})
    
    truncate=False
    if 'Lmax' in dictout:
        if 0 <= nmax < dictout['Lmax']:
            truncate = True

    #read side description data
    side1_d=np.fromfile(fid,dtype='|S24',count=nval1).astype('|U24')

    if dictout["type"] in ['BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN']:
        dictout["blockind"]=np.fromfile(fid,dtype=endianness+'I',count=nblocks)
        #also add coordinates in the full matrix
        coords=getBDcoords(dictout,trans)


    #possibly read second side description
    if dictout['type'] in ['BDFULLV0','BDFULLVN','FULLSQV0','FULLSQVN']:
        if vnum <= 2.2:
            side2_d=side1_d
        else:
            side2_d=np.fromfile(fid,dtype='|S24',count=nval1).astype('|U24')

    elif dictout["type"] == "FULL2DVN":
        side2_d=np.fromfile(fid,dtype='|S24',count=nval2).astype('|U24')


    # read vectors
    if nvec >0:
        vec=np.fromfile(fid,dtype=endianness+"d",count=nvec*nval1).reshape((nval1,nvec),'F')

    # read packed matrix
    pack=np.fromfile(fid,dtype=endianness+'d',count=pval1*pval2)
    
    if type(file_or_obj) == str:
        #expliclity close (otherwise it's assumed the caller will handle closure)
        fid.close()

    #possibly convert a SH side description in a (multi)-index
    try:
        mshi=get_shmi(side1_d)

        if truncate:
            # Truncate the input
            #create an index vector restricting the maximum degree of the output
            mshi2=xr.DataArray(mshi)
            #index corresponding to the rows and columns
            indx=mshi2.n <= nmax
            #index vector within the packed array
            pindx=indx[coords[0,:]].data*indx[coords[1,:]].data
            
            #get the subsets of the data
            side1_d=side1_d[indx]
            #new version
            mshi=get_shmi(side1_d)
            pack=pack[pindx]
            oldcoords=coords[:,pindx]
            idxnew=indx.cumsum()-1
            coords=np.zeros(oldcoords.shape,dtype=np.int64)
            coords[0,:]=idxnew[oldcoords[0,:]]
            coords[1,:]=idxnew[oldcoords[1,:]]
            nval1=len(mshi)

    except KeyError:
        pass
    
    if dictout["type"] in ['BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN']:
    #unpack in sparse matrix
        mat = sparse.COO(coords, pack, shape=(nval1,nval1),fill_value=0.0)
        # mat=sparse.GCXS.from_coo(mat,compressed_axes=[1])
    else: 
        raise NotImplemented(f"Cannot Unpack a matrix of {dictout['type']}")
    
    shname=SHindexBase.name
    shname_t=shname+"_"
    dsout=xr.Dataset(dict(mat=([shname,shname_t],mat)),coords=dict(nm=([shname],mshi),nm_=([shname_t],SHindexBase.mi_toggle(mshi))),attrs=dictout)
    #transform the data in dask 
    # dsout=dsout.chunk()
    if nvec > 0:
        dsout["vec"]=([shname],vec)

    return dsout

