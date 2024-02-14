# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++
# cython: profile=False

cimport cython
from gaunt cimport Gaunt,GauntReal

from cython.operator cimport dereference as deref
# import numpy as np 
cimport numpy as np
import xarray as xr
import pandas as pd

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
