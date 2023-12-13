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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Analysis:
    cdef public object _dsobj
    def __cinit__(self, lon=np.arange(-180.0,180.0,1.0), lat=np.arange(-90.0,90.0,1.0),grid=True):
        
        if type(lon) == xr.DataArray and type(lat) == xr.DataArray:
            #possibly check whether they share the same dimension name and size and set grid to false in that case
            if lon.dims[0] == lat.dims[0] and lon.size == lat.size:
                grid=False
            lon=lon.data
            lat=lat.data

        if grid:
            #create a xrray dataset which spans the in and output
            coords={"lat":("lat",lat),"lon":("lon",lon)}
            self._dsobj=xr.Dataset(coords=coords)
        else:
            if len(lon) != len(lat):
                raise RuntimeError("Analysis on a one dimensional set of points requires equally sized lon and lat ")

            coords={"lat":("nlonlat",lat),"lon":("nlonlat",lon)}
            self._dsobj=xr.Dataset(coords=coords)
                     
    def __call__(self,dain:xr.DataArray):
        """Perform the spherical harmonic analysis on an input xarray DataArray object""" 
        if type(dain) != xr.DataArray:
            raise RuntimeError("input type should be a xarray.DataArray")

        if "shi" not in dain.indexes:
            raise RuntimeError("Spherical harmonic index not found in input, cannot apply analysis operator to object")
        shidim= dain.get_axis_num('shi') 
        if shidim != dain.ndim-1 and shidim != 0:
            raise RuntimeError ("shi dimension must either be the first or last in input")
       
        # if not (dain.data.flags['C_CONTIGUOUS'] or dain.data.flags['F_CONTIGUOUS']):
            # raise RuntimeError("Cannot work with non-contiguous input arrays yet")

        coordsout={ky:val for ky,val in dain.coords.items() if val.ndim == 1  and val.dims[0] != "shi"}
        coordsout.update({ky:val for ky,val in self._dsobj.coords.items()})
        
        dimsin=[(dim,sz) for dim,sz in dain.sizes.items() if dim != "shi"]
        #it is import to have lat, and lon as first dimensions (slowest varying index)
        dimsout=[(dim,sz) for dim,sz in self._dsobj.sizes.items()]+dimsin

        #allocate space (c-contiguous) for the output data
        daout=xr.DataArray(np.zeros([val[1] for val in dimsout]),coords=coordsout,dims=[val[0] for val in dimsout])
        self._apply_dgemv(dain,daout)
        return daout

    # def _ddot_grid(self,dain:xr.DataArray,daout:xr.DataArray):
        # cdef int nmax=dain.sh.nmax
        # cdef int nmin=dain.sh.nmin
    
        # cdef double[::1] lonv=self._dsobj.lon.values
        # cdef double[::1] latv=self._dsobj.lat.values
        # nlat=len(latv)
        # nlon=len(lonv)

        # cdef Ynm ynm=Ynm(dain.indexes["shi"])
        
        

        # cdef double [:,::1] outv=daout.data
        # cdef double [::1] inval=dain.data
        
        # cdef:
            # int nsh=len(dain.indexes["shi"])
            # int incx=1
            # int incy=1
        
        # for ilat in range(nlat):
            # for ilon in range(nlon):
                # ynm.set(lonv[ilon],latv[ilat])
                # outv[ilat,ilon]=ddot(&nsh,&ynm.data[0],&incx,&inval[0],&incy)
        


    def _apply_dgemv(self,dain:xr.DataArray,daout:xr.DataArray):
        cdef int nmax=dain.sh.nmax
        cdef int nmin=dain.sh.nmin
        cdef double[::1] lonv=self._dsobj.lon.data.astype(np.double)
        cdef double[::1] latv=self._dsobj.lat.data.astype(np.double)
        cdef int nlat=len(latv)
        cdef int nlon=len(lonv)
        cdef bint grid

        
        cdef int shsize=len(dain.indexes['shi'])
        # product of the non-shi dimension lengths
        # note auxsize will be 1 (good!) when no other dimensions are present except for shi
        cdef int auxsize=np.prod([val for ky,val in dain.sizes.items() if ky != "shi"])
        
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
        
        elif  dain.get_axis_num('shi') == dain.ndim -1:
            #shi dimension is last
            inval=dain.data.reshape([auxsize,shsize]) # keep the sh dimension last
        elif  dain.get_axis_num('shi') == 0:
            #transpose input so the SHI dimension is last
            inval=dain.data.reshape([shsize,auxsize]).T # transpose to make the sh dimension last
        else:
            raise RuntimeError("Can only handle input where the shi dimension is first or last")
    
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


        # if (shIslast and dain.data.flags['C_CONTIGUOUS']) or (not shIslast and dain.data.flags['F_CONTIGUOUS']):
            # # Note this means that the shi dimension varies quickest
            # # For dgemv this means we need to transpose the input
            # # This can also be interpreted as a Fortran array where the rows span the shi dimension 
            # trans=b'T'
            # lda=shsize
            # m=shsize
            # n=auxsize

        # else:

            # #Note: shi dimension varies slowest
            # trans=b'N'
            # lda=auxsize
            # m=auxsize
            # n=shsize
        

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
        
        cdef Ynm ynm=Ynm(dain.indexes["shi"])
        cdef int ilat,ilon,idx          
        if grid:
            #use gridded approach (can be significantly faster because the latitudes will be in the outer loop)
            for ilat in range(nlat):
                for ilon in range(nlon):
                    ynm.set(lonv[ilon],latv[ilat])
                    dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,&ynm.data[0],&incx,&beta,&outv[ilat*nlon+ilon,0],&incy)

        else:
            for idx in range(npoints):
                ynm.set(lonv[idx],latv[idx])
                dgemv(&trans,&m,&n,&alpha,&inval[0,0],&lda,&ynm.data[0],&incx,&beta,&outv[idx,0],&incy)


