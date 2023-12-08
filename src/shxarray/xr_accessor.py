# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from .shxarbase import ShXrBase
from shxarray.earth.snrei import SnreiFactory

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
    
    def analysis(self,lon,lat,forcegrid=False,engine="shlib"):
        """
        Apply spherical harmonic analysis on a set of longitude, latitude points
        :param lon: Longitude in degrees East
        :param lat: Latitude in degrees North
        :param forcegrid: Force mapping to a grid when the longitude latitude have the same dimension/length
        :param engine: Spherical harmonic compute engine to use for the computation
        :return: A datarray for which the spherical harmonic coefficietn dimension is mapped to set of points
        The following scenarios can be handled:
        
        1: lon, lat are Xarray coordinate variables sharing the same di
        mension. Map to a list of points (SH dimension is mapped to a single dimension)
        2: lon, lat are Xarray coordinate variables with different dimensions: Map tot a grid spanned by lon,lat
        3 lon, lat are list-like objects with the same length: Map to a list of points unless (forcegrid=True)
        4. lon,lat are list-like objects of different lengths: Map to a grid
        """
        #dispatch to compute engine
        eng=self._eng(engine)()
        return eng.analysis(self._obj,lon,lat,forcegrid)
    

@xr.register_dataset_accessor("sh")
class SHDsAccessor(ShXrBase):
    snrei=SnreiFactory
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
