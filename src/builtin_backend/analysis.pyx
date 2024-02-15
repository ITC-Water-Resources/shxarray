# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

# distutils: language = c++
# cython: profile=True

import xarray as xr
import cython
# import numpy as np
cimport numpy as np
from cython.parallel cimport parallel, prange
from legendre cimport Ynm_cpp
from libc.stdio cimport printf
from openmp cimport omp_lock_t,omp_init_lock,omp_set_lock,omp_unset_lock
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.cf import find_lon,find_lat
from scipy.linalg.cython_blas cimport dger
from libc.math cimport cos

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Analysis:
    cdef public object _dsobj
    def __cinit__(self, int nmax):
        
        #create a spherical harmonic index
        self._dsobj=xr.Dataset(coords=SHindexBase.nm(nmax,0))
                     
    def __call__(self,dain:xr.DataArray):
        """Perform  spherical harmonic analysis on an input xarray DataArray object""" 
        if type(dain) != xr.DataArray:
            raise RuntimeError("input type should be a xarray.DataArray")

        #check for longitude/latitude grid and  grid separation
        loninfo=find_lon(dain.coords)  
        latinfo=find_lat(dain.coords)  
        if loninfo.step is None or latinfo.step is None:
            raise RuntimeError("grid is not equidistant in lon and lat directions")

        dain=dain.rename({loninfo.var.name:"lon",latinfo.var.name:"lat"})

        coordsout={ky:val for ky,val in dain.coords.items() if val.ndim == 1  and val.dims[0]  not in ["lon","lat"]}

        # add the SH index to the output coordinates
        coordsout.update({ky:val for ky,val in self._dsobj.coords.items()})
        
        dimsin=[(dim,sz) for dim,sz in dain.sizes.items() if dim not in ["lon","lat"]]
        # make sure the SHI index is the one who varies quickest (last dimension)
        dimsout=dimsin +[(dim,sz) for dim,sz in self._dsobj.sizes.items()]

        #allocate space (c-contiguous) for the output data
        daout=xr.DataArray(np.zeros([val[1] for val in dimsout]),coords=coordsout,dims=[val[0] for val in dimsout])
        
        #scale input with area weights
        cdef double stepx=loninfo.step*np.pi/180
        cdef double stepy=latinfo.step*np.pi/180
        # printf("stepx %f stepy %f\n",stepx,stepy)
        cdef double weight=stepx*stepy/(4*np.pi)
        self._apply_ana(dain,daout,weight)
        return daout
    
    cdef _apply_ana(self,dain:xr.DataArray,dout:xr.DataArray,double weight):


        cdef double[::1] lonv=dain.lon.data.astype(np.double)
        cdef double[::1] latv=dain.lat.data.astype(np.double)
        cdef int nlat=len(latv)
        cdef int nlon=len(lonv)

        cdef int auxsize=np.prod([val for ky,val in dain.sizes.items() if ky not in ["lon","lat"]])

        cdef int shsize=len(self._dsobj.indexes['nm'])
        #memoryview to output data (sh dimension should vary quickest) 
        cdef double [:,:] outv=dout.data.reshape([auxsize,shsize])
        #This is the same a s a Fortran contiguous array with dimension shsize,auxsize, lada=shsize  

        


        cdef int[::1] nv = self._dsobj.nm.n.data.astype(np.int32)
        cdef int[::1] mv = self._dsobj.nm.m.data.astype(np.int32)

        cdef Ynm_cpp[double] ynm
        cdef int ilat,ilon          

        cdef:
            int m=shsize
            int n=auxsize
            double alpha=1.0
            int lda=shsize
            int incx=1
            int incy=1
        
        cdef double [:,:,:] inval


        #check the order of lon lat 
        londim=dain.get_axis_num('lon') 
        latdim=dain.get_axis_num('lat')

        if abs(londim -latdim) != 1:
            raise RuntimeError("Longitude and latitude dimensions of input are not neighbouring, cannot handle this layout")
        
        cdef bint latfirst=False
        if londim == 0 or latdim == 0:
            if londim < latdim:
                inval=dain.data.reshape([nlon,nlat,auxsize],order='C')
                latfirst=False
            else:
                inval=dain.data.reshape([nlat,nlon,auxsize],order='C')
                latfirst=True
        elif londim == dain.ndim -1 or latdim == dain.ndim-1:
            if londim < latdim:
                inval=dain.data.reshape([auxsize,nlon,nlat],order='C').T
                latfirst=True
            else:
                inval=dain.data.reshape([auxsize,nlat,nlon],order='C').T
                latfirst=False
        cdef double d2r =np.pi/180
        incy=int(inval.strides[2]/8)
        cdef omp_lock_t lock
        omp_init_lock(&lock)
        with nogil, parallel():

            ynm=Ynm_cpp[double](shsize,&nv[0],&mv[0])
            for ilat in prange(nlat):
                alpha=weight*cos(latv[ilat]*d2r)
                for ilon in range(nlon):
                    ynm.set(lonv[ilon],latv[ilat])
                    # void dger(int *m, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *a, int *lda)   
                    #critical section for openmp because outv gets updated by all threads
                    omp_set_lock(&lock)
                    if latfirst:
                        dger(&m,&n,&alpha,ynm.data(),&incx,&inval[ilat,ilon,0],&incy,&outv[0,0],&lda)
                    else:
                        dger(&m,&n,&alpha,ynm.data(),&incx,&inval[ilon,ilat,0],&incy,&outv[0,0],&lda)
                    omp_unset_lock(&lock)

        



