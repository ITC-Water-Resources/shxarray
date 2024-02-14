# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++
# cython: profile=False

cimport cython
from wigner3j cimport Wigner3j

from cython.operator cimport dereference as deref
# import numpy as np 
cimport numpy as np
import xarray as xr
import pandas as pd

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def getWigner3j(j2,j3,m2,m3):
    """
    Compute non-zero Wigner3J symbols with their valid (j1,m1) for j2,j3,m2,m3 input
    """
    w3j = Wigner3j[double](j2,j3,m2,m3)
    assert(w3j.jmin() <= w3j.jmax())
    m=w3j.m()
    jm=pd.MultiIndex.from_tuples([(j,m) for j in range(w3j.jmin(),w3j.jmax()+1)],names=("j","m")) 
    return xr.DataArray(w3j.get(), coords=dict(jm=jm),dims=["jm"])

