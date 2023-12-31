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
from legendre cimport Ynm_cpp
from libc.stdio cimport printf
# from warnings import warn
# Todo: potentially improve speed by directly calling dgemv
# from scipy.linalg.blas import dgemv,ddot
from scipy.linalg.cython_blas cimport dgemv,ddot
from cython.operator cimport dereference as deref
from openmp cimport omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
from cpython.mem cimport PyMem_Malloc, PyMem_Free

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Analysis:
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
                
        # cdef Ynm ynm=Ynm(dain.indexes["shi"])
        cdef int[::1] nv = dain.shi.n.data.astype(np.int32)
        cdef int[::1] mv = dain.shi.m.data.astype(np.int32)
        cdef int[::1] tv = dain.shi.t.data.astype(np.int32)

        cdef Ynm_cpp[double] ynm
        cdef int ilat,ilon,idx          
        
        with nogil, parallel():
            ynm=Ynm_cpp[double](shsize,&nv[0],&mv[0],&tv[0])
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

@cython.profile(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def threadedynm():
    
    cdef int num_threads
    cdef int ithread
    cdef int max_threads=omp_get_max_threads()
    cdef np.ndarray ynms = np.empty([max_threads,],dtype=Ynm)
    cdef double * ynmdata
    cdef double [::1] lat=np.arange(-90.0,90.0,0.5)
    cdef double [::1] lon=np.arange(-180.0,180.0,0.5)
    cdef int nlat=len(lat)
    cdef int nlon=len(lon)
    cdef int ilat,ilon

    printf("max threads %d\n",max_threads)
    ynmlist = []
    # cdef Ynm_cpp[double] * ynm
    cdef Ynm_cpp[double]  ynm
    with nogil, parallel():
        num_threads = omp_get_num_threads()
        ithread = omp_get_thread_num()
        
        ynm = Ynm_cpp[double](300)
        # ynm = new Ynm_cpp[double](150)
    # printf("num_threads  %d,threadid %d assigned nmax %d\n",num_threads,ithread,1)
    # printf("thread id %d, ynmdata address %p, nmax %d\n",ithread,ynm.data(),ynm.nmax())
        for ilat in prange(nlat):
            # printf("thread id %d, doing latitude %f\n",ithread,lat[ilat])
            
            for ilon in range(nlon):
                # deref(ynm).set(lon[ilon],lat[ilat])
                ynm.set(lon[ilon],lat[ilat])
        # with gil:
            # # PyMem_Free(ynm)
            # del ynm
            

    return ynms
    # printf("memory address of ynm %p\n",ynm_ptr)

def checkynm():
    nmax=20
    da=xr.DataArray.sh.ones(20)
    shi=da.sh.truncate(12,3).shi
    ynm1=Ynm(shi)
    
    cdef int[::1] nv = np.array([n for n,_,_ in ynm1._shindex.values]).astype(np.int32)
    cdef int[::1] mv = np.array([m for _,m,_ in ynm1._shindex.values]).astype(np.int32)
    cdef int[::1] tv = np.array([t for _,_,t in ynm1._shindex.values]).astype(np.int32)
    
    cdef Ynm_cpp[double] ynm2= Ynm_cpp[double](len(ynm1),&nv[0],&mv[0],&tv[0])
    
    # cdef Ynm_cpp[double] ynm2= Ynm_cpp[double](nmax)
    cdef lon=40.3
    cdef lat=20.5
    cdef cython.size_t idx
    cdef int i,n,m,t
    # ynm1.set(lon,lat)
    ynm2.set(lon,lat)
    for it in ynm2.getmn():
        m=it.m
        n=it.n
        idx=it.i
        printf("%d n:%d m:%d\n",idx,n,m)

    # for i,(n,m,t) in enumerate(ynm1._shindex.values):
        # idx=ynm2.idx(n,m,t)
        # printf("%d n:%d m:%d t:%d %f %f\n",i,n,m,t,ynm1.data[i],ynm2[idx])

    



