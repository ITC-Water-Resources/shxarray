# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

# distutils: language = c++
# cython: profile=False

import cython
cimport numpy as np
import xarray as xr
from legendre cimport Legendre_nm
from libc.math cimport sin,cos,pi
from cython.operator cimport dereference as deref
from libcpp.pair cimport pair
from shxarray.sh_indexing import SHindexBase
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc
from libc.stdio cimport printf
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Ynm:
    cdef int nmax
    cdef double[::1] pnmcache
    cdef double[::1] data
    cdef long [::1] _idx
    cdef double latprev;
    cdef cython.size_t sz
    cdef bint resort
    # cdef Legendre_nm[double]*legnm_ptr
    cdef Legendre_nm[double] legnm
    cdef public object _shindex
    def __cinit__(self,nmax_or_index):
        if type(nmax_or_index) == int:
            self.nmax=nmax_or_index
            self.sz=SHindexBase.nsh(self.nmax,squeeze=True)
            self.resort=False
        else:
            self._shindex=nmax_or_index
            self.nmax=nmax_or_index.max()[0]
            self.sz=len(nmax_or_index)
            self.resort=True


        cdef int[:,::1] nmt=np.zeros([self.sz,3],dtype=np.int32)

        cdef cython.size_t i=0
        for m in range(self.nmax+1):
            for n in range(m,self.nmax+1):
                nmt[i,0]=n
                nmt[i,1]=m
                nmt[i,2]=0
                i+=1
                if m > 0:
                    nmt[i,0]=n
                    nmt[i,1]=m
                    nmt[i,2]=1
                    i+=1


        
        if self.resort:
            # get the resorting index
            self._idx=self._shindex.get_indexer(SHindexBase.mi_fromarrays(np.asarray(nmt).T))
        else:
            #create a new index from the nmt array
            self._shindex=SHindexBase.mi_fromarrays(np.asarray(nmt).T)

        # self.legnm_ptr=new Legendre_nm[double](self.nmax)
        
        self.legnm=Legendre_nm[double](self.nmax)
        self.latprev=-1000.0 #an invalid initial value triggers initial computation
        
        # cdef cython.size_t pnm_sz=deref(self.legnm_ptr).size()
        cdef cython.size_t pnm_sz=self.legnm.size()
        #initialize cache and lookup index
        self.pnmcache=np.zeros([pnm_sz])
        self.data=np.zeros([self.sz])


    def __dealloc__(self):
        pass
        # del self.legnm_ptr

    @property
    def nmax(self):
        return self.nmax
    
    def __len__(self):
        return self.sz
    
    def  __call__(self,lon, lat):
        
        cdef int npos
        cdef double[:,::1] data;

        if np.isscalar(lon) and np.isscalar(lat):
            self.set(lon,lat)
            dsout=xr.DataArray(self.data,coords={"shi":self._shindex,"lon":lon,"lat":lat},dims=["shi"],name="Ynm")
        else:
            #multiple sets requested
            if len(lon) != len(lat):
                raise RuntimeError("input longitude and latitude needs to be of the same length")
            npos=len(lon)
            data=np.empty([npos,self.sz])
            
            for i in range(npos):
                self.set(lon[i],lat[i])
                data[i,:]=self.data
            
            dsout=xr.DataArray(data,coords={"shi":("shi",self._shindex),"lon":("nlonlat",lon),"lat":("nlonlat",lat)},dims=["nlonlat","shi"],name="Ynm")
        
        return dsout

    cdef void set(self,double lon,double lat) noexcept nogil:
        cdef double costheta
        if (lat != self.latprev):
            #(re)compute associated Legendre functions
            costheta=sin(lat*pi/180.0)
            self.latprev=lat
            # deref(self.legnm_ptr).set(costheta,&self.pnmcache[0])
            self.legnm.set(costheta,&self.pnmcache[0])
        cdef double lonr=lon*pi/180.0

        cdef int m=0
        cdef int n=0
        cdef double c_mlambda=0
        cdef double s_mlambda=0

        #note: the order of these loops and their corresponding n,m,t MUST agree with the index self._idx!
        cdef cython.size_t i=0
        cdef cython.size_t ipnm=0
        cdef long idx=0

        for m in range(self.nmax+1):
            c_mlambda=cos(m*lonr)
            s_mlambda=sin(m*lonr)
            for n in range(m,self.nmax+1):
                # ipnm=deref(self.legnm_ptr).idx(n,m)
                ipnm=self.legnm.idx(n,m)
                if self.resort:
                    idx=self._idx[i]
                else:
                    idx=i

                self.data[idx]=c_mlambda* self.pnmcache[ipnm]
                i+=1
                if m > 0:
                    if self.resort:
                        idx=self._idx[i]
                    else:
                        idx=i
                    self.data[idx]=s_mlambda* self.pnmcache[ipnm]
                    i+=1

                
