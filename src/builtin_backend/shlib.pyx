# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

include "legendre.pyx" 
include "wigner3j.pyx"
include "ynm.pyx"
include "analysis.pyx"

from shxarray.shcomputebase import SHComputeBackendBase
# from shxarray.shlib import Analysis

class SHComputeBackend(SHComputeBackendBase):
    def analysis(self,dain,lon,lat):
        ana=Analysis(lon,lat)
        return ana(dain)




