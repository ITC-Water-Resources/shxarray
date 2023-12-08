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
from shxarray.shlib import Ynm
from libc.stdio cimport printf
# from warnings import warn
# Todo: potentially improve speed by directly calling dgemv
# from scipy.linalg.blas import dgemv,ddot
from scipy.linalg.cython_blas cimport dgemv,ddot
from cython.operator cimport dereference as deref

@cython.boundscheck(True)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Analysis:
    cdef public object _dsobj
    def __cinit__(self, lon=np.arange(-180.0,180.0,1.0), lat=np.arange(-90.0,90.0,1.0)):

        #create a xrray dataset which spans the in and output
        coords={"lat":("lat",lat),"lon":("lon",lon)}
        self._dsobj=xr.Dataset(coords=coords)
         
    def __call__(self,dain:xr.DataArray):
        """Perform the spherical harmonic analysis on an input xarray DataArray object""" 
        if type(dain) != xr.DataArray:
            raise RuntimeError("input type should be a xarray.DataArray")

        if "shi" not in dain.indexes:
            raise RuntimeError("Spherical harmonic index not found in input, cannot apply analysis operator to object")
        cdef int nmax=dain.sh.nmax
        cdef int nmin=dain.sh.nmin
        cdef double[::1] lonv=self._dsobj.lon.values
        cdef double[::1] latv=self._dsobj.lat.values
        nlat=len(latv)
        nlon=len(lonv)
        newcoords={ky:val for ky,val in dain.coords.items() if val.dims[0] != "shi"}
        newcoords.update({ky:val for ky,val in self._dsobj.coords.items()})
    
        dims=[(dim,sz) for dim,sz in self._dsobj.dims.items()]

        for dim,sz in dain.sizes.items():
            if dim != "shi":
                dims.append((dim,sz))

        #allocate space for the output data
        daout=xr.DataArray(np.zeros([val[1] for val in dims]),coords=newcoords,dims=[val[0] for val in dims])


        cdef double [:,::1] outv=daout.data
        
        cdef double [::1] inval=dain.data
        
        cdef:
            double alpha=1.0
            double beta=0.0
            int nsh=len(dain.indexes["shi"])
            int offx=0
            int offy=0
            int incx=1
            int incy=1
            int transa=0
        cdef Ynm ynm=Ynm(dain.indexes["shi"])

        # with nogil, parallel():
            # #make sure to generate surface spherical harmonics in the same order as the input 
            # ynm=Ynm(dain.indexes["shi"])
        for ilat in range(nlat):
            for ilon in range(nlon):
                ynm.set(lonv[ilon],latv[ilat])
                outv[ilat,ilon]=ddot(&nsh,&ynm.data[0],&incx,&inval[0],&incy)
        return daout
