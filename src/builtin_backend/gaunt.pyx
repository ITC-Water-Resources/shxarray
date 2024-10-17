# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++
# cython: profile=False

cimport cython
from gaunt cimport Gaunt,GauntReal

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map
cimport numpy as np
import xarray as xr
import pandas as pd
from libc.math cimport sqrt

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
    nm=pd.MultiIndex.from_tuples(gauntreal.nm(),names=("n","m")) 
    return xr.DataArray(gauntreal.get(), coords=dict(nm=nm),dims=["nm"])


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
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
    
    
     
    
    cdef vector[pair[int,int]] nm1
    cdef GauntReal[double] gauntreal
    
    cdef int n1,m1,n2,n3,m2,m3
    cdef cython.size_t ngnt,ix,iy,iz
    cdef double val
    
    #look up map for input coefficient
    cdef map[pair[int,int],cython.size_t] nm1_map
    cdef cython.size_t idx=0

    for (n1,m1) in daobj.nm.data:
        nm1_map[(n1,m1)]=idx
        idx+=1
    
    ngnt=0
    val=0.0

    with nogil, parallel():
        for ix in prange(nsh,schedule="guided"):
            n2=n[ix]
            m2=m[ix]
            for iy in range(0,ix+1):
                #only compute for the upper triangle
                n3=n[iy]
                m3=m[iy]
                gauntreal = GauntReal[double](n2,n3,m2,m3)
                ngnt=gauntreal.size()
                nm1=gauntreal.nm()
                val=0.0
                for iz in range(ngnt):
                    val=val+gauntreal[iz]*davec[nm1_map[nm1[iz]]]
                p2smat[ix,iy]=val*norm
                p2smat[iy,ix]=val*norm

    return dsp2s
