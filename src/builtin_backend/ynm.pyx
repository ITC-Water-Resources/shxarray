# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

# distutils: language = c++
# cython: profile=False

import cython
cimport numpy as np
import xarray as xr
from legendre cimport Ynm_cpp,mni
from shxarray.sh_indexing import SHindexBase
from libc.stdio cimport printf


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Ynm:
    cdef Ynm_cpp[double] _ynm
    cdef double[::1] data
    cdef public object _shindex
    def __cinit__(self,nmax_or_index):
        
        cdef int[::1] nv  
        cdef int[::1] mv 
        cdef int[::1] tv
        cdef int [:,::1] nmt
        cdef cython.size_t sz,idx
        cdef int nmax,n,m,t
        cdef mni it
        if type(nmax_or_index) == int:
            nmax=nmax_or_index
            self._ynm=Ynm_cpp[double](nmax)
            sz=self._ynm.size()
            #create a sh index
            nmt=np.zeros([sz,3],dtype=np.int32)
            for it in self._ynm.getmn():
                n=it.n
                m=it.m
                idx=it.i
                if m<0:
                    nmt[idx,0]=n
                    nmt[idx,1]=-m
                    nmt[idx,2]=1
                else:
                    nmt[idx,0]=n
                    nmt[idx,1]=m
                    nmt[idx,2]=0

            self._shindex=SHindexBase.mi_fromarrays(np.asarray(nmt).T)
        else:
            nv=np.array([n for n,_,_ in nmax_or_index.values]).astype(np.int32)
            mv=np.array([m for _,m,_ in nmax_or_index.values]).astype(np.int32)
            tv=np.array([t for _,_,t in nmax_or_index.values]).astype(np.int32)

            sz=len(nmax_or_index)
            self._ynm=Ynm_cpp[double](sz,&nv[0],&mv[0],&tv[0])
            self._shindex=nmax_or_index

        #have data memory view point to the memory of the cpp class
        self.data = <double[:sz:1]>(self._ynm.data()) 


    @property
    def nmax(self):
        return self._ynm.nmax()
    
    def __len__(self):
        return self._ynm.size()
    
    def  __call__(self,lon, lat):
        
        cdef int npos
        cdef double[:,::1] data;

        if np.isscalar(lon) and np.isscalar(lat):
            
            self._ynm.set(lon,lat)
            dsout=xr.DataArray(self.data,coords={"shi":self._shindex,"lon":lon,"lat":lat},dims=["shi"],name="Ynm")
        else:
            #multiple sets requested
            if len(lon) != len(lat):
                raise RuntimeError("input longitude and latitude needs to be of the same length")
            npos=len(lon)
            data=np.empty([npos,self._ynm.size()])
            
            for i in range(npos):
                self._ynm.set(lon[i],lat[i])
                data[i,:]=self.data
            
            dsout=xr.DataArray(data,coords={"shi":("shi",self._shindex),"lon":("nlonlat",lon),"lat":("nlonlat",lat)},dims=["nlonlat","shi"],name="Ynm")
        
        return dsout

