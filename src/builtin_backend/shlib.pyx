# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
"""
**shlib** is shxarray's default binary Cython backend. 
Some of the heavy lifting such as synthesis and analysis operations, is done using this the functions of this shared library.
"""

include "legendre.pyx" 
include "wigner3j.pyx"
include "gaunt.pyx"
include "ynm.pyx"
include "synthesis.pyx"
include "analysis.pyx"

from shxarray.core.shcomputebase import SHComputeBackendBase

class SHComputeBackend(SHComputeBackendBase):
    def synthesis(self,dain,lon, lat,grid):
        syn=Synthesis(lon,lat,grid)
        return syn(dain)


    def analysis(self,dain,nmax,method):
        if method not in ['integrate']:
            raise RuntimeError (f"SHComputeBackend does not support {method} in analysis")
        ana=Analysis(nmax)
        return ana(dain)
    
    def gaunt(self,n2,n3,m2,m3):
        return getGaunt(n2,n3,m2,m3)

    def gauntReal(self,n2,n3,m2,m3):
        return getGauntReal(n2,n3,m2,m3)

    def wigner3j(self,j2,j3,m2,m3):
        return getWigner3j(j2,j3,m2,m3)
