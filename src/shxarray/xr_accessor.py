# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from .shxarbase import ShXrBase



@xr.register_dataarray_accessor("sh")
class SHDaAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    
    @staticmethod
    def zeros(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={}):
        """0-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,squeeze,name,auxcoords)
    
    @staticmethod
    def ones(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={}):
        """1-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,squeeze,name,auxcoords)

@xr.register_dataset_accessor("sh")
class SHDsAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
        
    @staticmethod
    def zeros(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={}):
        """0-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,squeeze,name,auxcoords).to_dataset()
    
    @staticmethod
    def ones(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={}):
        """1-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,squeeze,name,auxcoords).to_dataset()
