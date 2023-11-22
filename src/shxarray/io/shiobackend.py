# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023

from xarray.backends import BackendEntrypoint
from shxarray.logging import logger 
from shxarray.io.icgem import readIcgem
from shxarray.io.gsmv6 import readGSMv6
import os
import re

class ICGEMBackEntryPoint(BackendEntrypoint):
    url="https://github.com/ITC-Water-Resources/shxarray"
    description = "Read spherical harmonic coefficients in ICGEM format"
    def open_dataset(self,filename_or_obj,*,drop_variables=None):
        dsout=readIcgem(filename_or_obj)
        if drop_variables is not None:
            dsout=dsout.drop_vars(drop_variables)
        return dsout
    
    def guess_can_open(self,filename_or_obj):
        try:
            strrep=str(filename_or_obj)
            if strrep.endswith(".gfc"):
                return True
            if strrep.endswith("gfc.gz"):
                return True
        except AttributeError:
            return False
            
        return False


class GSMv6BackEntryPoint(BackendEntrypoint):
    url="https://github.com/ITC-Water-Resources/shxarray"
    description = "Read spherical harmonic coefficients in GSM-V6 format"
    def open_dataset(self,filename_or_obj,*,drop_variables=None):
        dsout=readGSMv6(filename_or_obj)
        if drop_variables is not None:
            dsout=dsout.drop_vars(drop_variables)
        return dsout
    
    def guess_can_open(self,filename_or_obj):
        try:
            strrep=str(filename_or_obj)
            # search for conventional file naming of GRACE
            if re.search('G[SA][MBCA]-2[^\s]*.gz',strrep):
                return True
        except AttributeError:
            return False
            
        return False


