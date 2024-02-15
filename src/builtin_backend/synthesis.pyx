# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

# distutils: language = c++
# cython: profile=False

import xarray as xr
import cython
cimport numpy as np
from cython.parallel cimport parallel, prange
from legendre cimport Ynm_cpp
from libc.stdio cimport printf
# from warnings import warn
from scipy.linalg.cython_blas cimport dgemv
from shxarray.core.cf import get_cfatts
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Synthesis:
    cdef public object _dsobj
    def __cinit__(self, lon, lat,grid=True):
        
        if type(lon) == xr.DataArray and type(lat) == xr.DataArray:
            #possibly check whether they share the same dimension name and size and set grid to false in that case
            if lon.dims[0] == lat.dims[0] and lon.size == lat.size:
                grid=False
            lon=lon.data
            lat=lat.data

        if grid:
            #create a xrray dataset which spans the in and output
            coords={"lat":("lat",lat,get_cfatts('latitude')),"lon":("lon",lon,get_cfatts('longitude'))}
            self._dsobj=xr.Dataset(coords=coords)
        else:
            if len(lon) != len(lat):
                raise RuntimeError("Synthesis on a one dimensional set of points requires equally sized lon and lat ")

            coords={"lat":("nlonlat",lat,get_cfatts('latitude')),"lon":("nlonlat",lon,get_cfatts('longitude'))}
            self._dsobj=xr.Dataset(coords=coords)

                     
    def __call__(self,dain:xr.DataArray):
        """Perform the spherical harmonic synthesis on an input xarray DataArray object""" 
        if type(dain) != xr.DataArray:
            raise RuntimeError("input type should be a xarray.DataArray")

        if "nm" not in dain.indexes:
            raise RuntimeError("Spherical harmonic index not found in input, cannot apply synthesis operator to object")
        shidim= dain.get_axis_num('nm') 
        if shidim != dain.ndim-1 and shidim != 0:
            raise RuntimeError ("nm dimension must either be the first or last in input")
       
        # if not (dain.data.flags['C_CONTIGUOUS'] or dain.data.flags['F_CONTIGUOUS']):
            # raise RuntimeError("Cannot work with non-contiguous input arrays yet")

        coordsout={ky:val for ky,val in dain.coords.items() if val.ndim == 1  and val.dims[0] != "nm"}
        coordsout.update({ky:val for ky,val in self._dsobj.coords.items()})
        
        dimsin=[(dim,sz) for dim,sz in dain.sizes.items() if dim != "nm"]
        #it is import to have lat, and lon as first dimensions (slowest varying index)
        dimsout=[(dim,sz) for dim,sz in self._dsobj.sizes.items()]+dimsin

        #allocate space (c-contiguous) for the output data
        daout=xr.DataArray(np.zeros([val[1] for val in dimsout]),coords=coordsout,dims=[val[0] for val in dimsout])
        self._apply_dgemv(dain,daout)
        return daout


    cdef _apply_dgemv(self,dain:xr.DataArray,daout:xr.DataArray):
        cdef int nmax=dain.sh.nmax
        cdef int nmin=dain.sh.nmin
        cdef double[::1] lonv=self._dsobj.lon.data.astype(np.double)
        cdef double[::1] latv=self._dsobj.lat.data.astype(np.double)
        cdef int nlat=len(latv)
        cdef int nlon=len(lonv)
        cdef bint grid

        
        cdef int shsize=len(dain.indexes['nm'])
        # product of the non-shi dimension lengths
        # note auxsize will be 1 (good!) when no other dimensions are present except for shi
        cdef int auxsize=np.prod([val for ky,val in dain.sizes.items() if ky != "nm"])
        
        cdef int npoints
        if self._dsobj.coords['lat'].dims[0] == self._dsobj.coords['lon'].dims[0]:
            #note lon and lat share the same dimensions
            npoints=len(self._dsobj.coords['lat'])
            grid=False
        else:
            grid=True
            npoints=nlon*nlat

        
        # printf("auxsize,shsize %d %d\n",auxsize,shsize)

        
        cdef double [:,:] inval
        if dain.ndim == 1:
            if dain.data.strides[0]/8 == 1:
                inval=dain.data.reshape([auxsize,shsize]) 
            else:
                inval=dain.data.reshape([shsize,auxsize]).T 
        
        elif  dain.get_axis_num('nm') == dain.ndim -1:
            #nm dimension is last
            inval=dain.data.reshape([auxsize,shsize]) # keep the sh dimension last
        elif  dain.get_axis_num('nm') == 0:
            #transpose input so the SHI dimension is last
            inval=dain.data.reshape([shsize,auxsize]).T # transpose to make the sh dimension last
        else:
            raise RuntimeError("Can only handle input where the nm dimension is first or last")
    
        # cdef int sz1,sz2
        # sz1=int(inval.strides[0]/8)
        # sz2=int(inval.strides[1]/8)
        # sz1=dain.data.strides[0]/8
        # sz2=dain.data.strides[1]/8

        # printf("Strides %d %d\n",sz1,sz2)
        


        cdef:
            char trans=b'N'
            int m
            int n
            double alpha=1.0
            double beta=0.0
            int lda
            int incx=1
            int incy=1

        cdef double [:,:] outv=daout.data.reshape([npoints,auxsize])
        
        if inval.strides[1] < inval.strides[0]:
            # printf("Transpose call\n")
            trans=b'T'
            lda=shsize
            m=shsize
            n=auxsize
        else:
            # printf("Normal call\n")
            trans=b'N'
            lda=int(inval.strides[1]/8)
            m=auxsize
            n=shsize


        #fortran call signature dgemv( 	character  	trans, 
        #integer  	m,
        #integer  	n,
        #double precision  	alpha,
        #double precision, dimension(lda,*)  	a,
        #integer  	lda,
        #double precision, dimension(*)  	x,
        #integer  	incx,
        #double precision  	beta,
        #double precision, dimension(*)  	y,
        #integer  	incy 
                
        cdef int[::1] nv = dain.nm.n.data.astype(np.int32)
        cdef int[::1] mv = dain.nm.m.data.astype(np.int32)
        # cdef int[::1] tv = dain.shi.t.data.astype(np.int32)

        cdef Ynm_cpp[double] ynm
        cdef int ilat,ilon,idx          
        
        with nogil, parallel():
            ynm=Ynm_cpp[double](shsize,&nv[0],&mv[0])
            if grid:
                #use gridded approach (can be significantly faster because the latitudes will be in the outer loop)
                for ilat in prange(nlat):
                    for ilon in range(nlon):
                        ynm.set(lonv[ilon],latv[ilat])
                        # dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,&ynm.data[0],&incx,&beta,&outv[ilat*nlon+ilon,0],&incy)
                        dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,ynm.data(),&incx,&beta,&outv[ilat*nlon+ilon,0],&incy)

            else:
                for idx in prange(npoints):
                    ynm.set(lonv[idx],latv[idx])
                    # dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,&ynm.data[0],&incx,&beta,&outv[idx,0],&incy)
                    dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,ynm.data(),&incx,&beta,&outv[idx,0],&incy)

