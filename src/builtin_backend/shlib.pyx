# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

include "legendre.pyx" 
include "wigner3j.pyx"
include "ynm.pyx"
include "synthesis.pyx"
include "analysis.pyx"

from shxarray.shcomputebase import SHComputeBackendBase

class SHComputeBackend(SHComputeBackendBase):
    def synthesis(self,dain,lon, lat,grid):
        syn=Synthesis(lon,lat,grid)
        return syn(dain)


    def analysis(self,dain,nmax,method):
        if method not in ['integrate']:
            raise RuntimeError (f"SHComputeBackend does not support {method} in analysis")
        ana=Analysis(nmax)
        return ana(dain)


