# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++

cimport cython
from gaunt cimport Gaunt,GauntReal
from shindex cimport Nmindex

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map
cimport numpy as np
import xarray as xr
import pandas as pd
from libc.math cimport sqrt
from openmp cimport omp_lock_t,omp_init_lock,omp_set_lock,omp_unset_lock,omp_get_max_threads,omp_get_thread_num

from cython.operator cimport dereference as deref

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def getGaunt(n2,n3,m2,m3):
    """
    Compute non-zero Gaunt coefficients for valid values of n1 and m1
    """
    gaunt = Gaunt[double](n2,n3,m2,m3)
    assert(gaunt.nmin() <= gaunt.nmax())
    m=gaunt.m()
    nm=pd.MultiIndex.from_tuples([(n,m) for n in range(gaunt.nmin(),gaunt.nmax()+1,2)],names=("n","m")) 
    return xr.DataArray(gaunt.get(), coords=dict(nm=nm),dims=["nm"])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def getGauntReal(n2,n3,m2,m3):
    """
    Compute non-zero Real Gaunt coefficients for valid values of n1 and m1
    """
    gauntreal = GauntReal[double](n2,n3,m2,m3)
    assert(gauntreal.nmin() <= gauntreal.nmax())
    nm=pd.MultiIndex.from_tuples(gauntreal.nmvec(),names=("n","m")) 
    return xr.DataArray(gauntreal.get(), coords=dict(nm=nm),dims=["nm"])


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.profile(False)
def getp2s(daobj):
    #Initialize the output array with zeros
    cdef int nmax=int(daobj.sh.nmax/2)
    dsp2s=xr.DataArray.sh.zeros(nmax,nshdims=2)
    cdef int nsh=len(dsp2s['nm'])
    cdef double norm=sqrt(4*np.pi)
    
    cdef double[:,:] p2smat=dsp2s.data#.astype(np.double)
    
    cdef double[:] davec=daobj.data.astype(np.double)

    cdef int[::1] n = dsp2s.nm.n.data.astype(np.int32)
    cdef int [::1] m = dsp2s.nm.m.data.astype(np.int32)
    
    
     
    
    # cdef vector[pair[int,int]] nm1
    cdef GauntReal[double] gauntreal
    cdef int n1,m1,n2,n3,m2,m3
    cdef cython.size_t ngnt,ix,iy,iz
    cdef double val
    
    #quick look up map for input coefficient
    #cdef map[pair[int,int],cython.size_t] nm1_map
    cdef Nmindex nm1_map=Nmindex(len(daobj))
    cdef cython.size_t idx=0

    for (n1,m1) in daobj.nm.data:
        # nm1_map[(n1,m1)]=idx
        #nm1_map[(n1,m1)]=idx
        nm1_map.set((n1,m1),idx)
        idx+=1


    
    ngnt=0
    val=0.0

    with nogil, parallel():
        gauntreal = GauntReal[double](nmax*2)
        for ix in prange(nsh,schedule="guided"):
            n2=n[ix]
            m2=m[ix]
            for iy in range(0,ix+1):
                #only compute for the upper triangle
                n3=n[iy]
                m3=m[iy]
                #gauntreal = GauntReal[double](n2,n3,m2,m3)
                gauntreal.set(n2,n3,m2,m3)
                ngnt=gauntreal.size()
                val=0.0
                for iz in range(ngnt):
                    val=val+gauntreal[iz]*davec[nm1_map[gauntreal.nm(iz)]]
                p2smat[ix,iy]=val*norm
                p2smat[iy,ix]=val*norm

    return dsp2s


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.profile(False)
def multiply_sh(dash1,dash2):
    """
    Compute the multiplication of two spherical harmonic datasets while staying in the spectral domain  
    """

    cdef int shidim1=dash1.get_axis_num('nm')

    cdef int shidim2=dash2.get_axis_num('nm')

    if shidim1 != dash1.ndim - 1 or shidim2 != dash2.ndim - 1:
        raise ValueError("Input data should have the spherical harmonic dimension as the last dimension, consider transposing the input")



    cdef int nmax1=dash1.sh.nmax
    cdef int nmax2=dash2.sh.nmax
    cdef int nmaxout=nmax1+nmax2

    # if nmax1 < nmax2:
        # raise ValueError("nmax of left input system should be larger or equal than the second input system, consider swopping the input")

    
    cdef cython.size_t auxsize1,auxsize2
    # cdef cython.size_t auxsize1= np.prod([val for ky,val in dash1.sizes.items() if ky not in ["nm"]])
    
    # cdef cython.size_t auxsize2= np.prod([val for ky,val in dash2.sizes.items() if ky not in ["nm"]])
    
    auxdims=[dim for dim in dash1.dims if dim not in ["nm"]]

    if len(auxdims) > 1:
        dash1=dash1.stack(auxdim=auxdims)
        dash2=dash2.stack(auxdim=auxdims)
        auxdim='auxdim_'
        
    elif len(auxdims) == 0:
        auxdim="auxsingle_"
        #expand for consistency with the below code
        dash1=dash1.expand_dims(auxdim)
        dash2=dash2.expand_dims(auxdim)
    else:
        auxdim=auxdims[0]
    
    auxsize1=dash1.sizes[auxdim]    
    auxsize2=dash2.sizes[auxdim]    
    

    if auxsize1 != auxsize2:
        raise ValueError("Auxiliary dimension size should be the same for both input systems")
    


    outcoords=SHindexBase.nm(nmaxout,0)
    cdef int nsh3=SHindexBase.nsh(nmaxout,0)
    auxcoords1=dash1.sh.auxcoords()

    if auxcoords1 and auxdim in auxcoords1:
        outcoords[auxdim]=auxcoords1[auxdim]

    dimsout=[(auxdim,auxsize1),("nm",nsh3)]
    
    #Initialize the output array with zeros
    #make sure that the sh dimension is last

    daout=xr.DataArray(np.zeros([vl[1] for vl in dimsout],dtype=np.double,order="C"),dims=[vl[0] for vl in dimsout],coords=outcoords)


    cdef int nsh1=dash1.sizes['nm']
    cdef int nsh2=dash2.sizes['nm']


    cdef double norm=sqrt(4*np.pi)
    
    cdef double[:,:] dash1mat

    if dash1.data.strides[0] == 8:
        #dash1 is C-contiguous
        dash1mat=dash1.data.reshape([auxsize1,nsh1])
    else:
        #dash1 is Fortran contiguous
        dash1mat=dash1.data.reshape([nsh1,auxsize1]).T
    
    cdef double[:,:] dash2mat

    if dash2.data.strides[0] == 8:
        #dash1 is C- contiguous
        dash2mat=dash2.data.reshape([auxsize2,nsh2])
    else:
        #dash1 is Fortran contiguous
        dash2mat=dash2.data.reshape([nsh2,auxsize2]).T


    cdef int[::1] nd1 = dash1.nm.n.data.astype(np.int32)
    cdef int [::1] md1 = dash1.nm.m.data.astype(np.int32)
    
    cdef int[::1] nd2 = dash2.nm.n.data.astype(np.int32)
    cdef int [::1] md2 = dash2.nm.m.data.astype(np.int32)
    
     
    
    cdef GauntReal[double] gauntreal
    cdef int n1,m1,n2,n3,m2,m3
    cdef cython.size_t ix,iy,iz,i1t,i2t,iaux,i2
    cdef double cccontrib
    
    #quick look up map for the output coefficient
    cdef cython.size_t ngnt=0
    cdef Nmindex nm3_map=Nmindex(nsh3)
    cdef pair[int,int] nm3 
    cdef cython.size_t idx=0
    nm3_map=Nmindex(nsh3)

    for (n3,m3) in daout.nm.data:
        nm3_map.set(pair[int,int](n3,m3),idx)
        idx+=1
    cdef int nthreads=omp_get_max_threads()
    cdef double[:,:,:] daoutcompound=np.zeros([nthreads,auxsize1,nsh3],dtype=np.double,order="C")
    cdef int ithread  
    with nogil, parallel():
        ithread=omp_get_thread_num()
        # with gil:
            # shxlogger.info("Thread %d of %d started" % (ithread,nthreads))
        gauntreal = GauntReal[double](nmaxout)
        for ix in prange(nsh1,schedule="guided"):
        # for ix in range(nsh1):
            n1=nd1[ix]
            m1=md1[ix]
            for iy in range(nsh2):
                #only needed to compute Real Gaunt coefficients for one combination due to symmetry
                n2=nd2[iy]
                m2=md2[iy]

                gauntreal.set(n1,n2,m1,m2)
                ngnt=gauntreal.size()
                for iaux in range(auxsize1):
                    cccontrib=dash1mat[iaux,ix]*dash2mat[iaux,iy]
                    for iz in range(ngnt):
                        nm3=gauntreal.nm(iz)
                        idx=nm3_map[gauntreal.nm(iz)]
                        daoutcompound[ithread,iaux,idx]+=gauntreal[iz]*(cccontrib)*norm

    #Sum the results from all threads (note the [:,:] syntax to ensure an inplace copy)   
    daout[:,:]=np.sum(daoutcompound,axis=0)
    #unstack the auxiliary dimensions if needed 
    if auxdim == 'auxdim_':
        daout=daout.unstack(auxdim)
    elif auxdim == 'auxsingle_':
        daout=daout.squeeze(auxdim,drop=True)
    # omp_destroy_lock(&lock)
    return daout




