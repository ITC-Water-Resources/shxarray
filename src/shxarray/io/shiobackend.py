# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023

from xarray.backends import BackendEntrypoint
from shxarray.core.logging import logger 
from shxarray.io.icgem import readIcgem
from shxarray.io.gsmv6 import readGSMv6
from shxarray.io.binv_legacy import readBINV
from shxarray.io.shascii import readSHAscii
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
            if re.search(r'G[SA][MBCA]-2[^\s]*.gz',strrep):
                return True
        except AttributeError:
            return False
            
        return False

class SHAsciiBackEntryPoint(BackendEntrypoint):
    url="https://github.com/ITC-Water-Resources/shxarray"
    description = "Read spherical harmonic coefficients in generic n,m, cnm, snm, sigcnm, sigsnm ascii format"
    def open_dataset(self,filename_or_obj,*,drop_variables=None):
        dsout=readSHAscii(filename_or_obj)
        if drop_variables is not None:
            dsout=dsout.drop_vars(drop_variables)
        return dsout
    
    def guess_can_open(self,filename_or_obj):
        #User need to use this engine explicitly as the filenaming can be anything
        return False

## NOte: this currently does not work (xarray changes the underlying sparse array)
class DDKBackEntryPoint(BackendEntrypoint):
    url="https://github.com/ITC-Water-Resources/shxarray"
    description = "Read spherical harmonic filter coefficients in legacy BINV format"
    def open_dataset(self,filename_or_obj,*,drop_variables=None):
        dsout=readBINV(filename_or_obj)
        breakpoint()
        if drop_variables is not None:
            dsout=dsout.drop_vars(drop_variables)
        return dsout
    
    def guess_can_open(self,filename_or_obj):
        try:
            strrep=str(filename_or_obj)
            # search for conventional file naming of  DDK files
            if strrep.startswith('Wbd_2-120'):
                #known anisotropic DDK filter matrix
                return True
        except AttributeError:
            return False
            
        return False
