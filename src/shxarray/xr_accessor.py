# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from .shxarbase import ShXrBase
from shxarray.earth.snrei import SnreiFactory
import numpy as np

@xr.register_dataarray_accessor("sh")
class SHDaAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    
    @staticmethod
    def zeros(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C'):
        """0-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,squeeze,name,auxcoords,order=order)
    
    @staticmethod
    def ones(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C'):
        """1-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,squeeze,name,auxcoords,order=order)
    
    def analysis(self,lon=np.arange(-180.0,180.0,1.0), lat=np.arange(-90.0,90.0,1.0),grid=True,engine="shlib"):
        """
        Apply spherical harmonic analysis on a set of longitude, latitude points
        :param lon: Longitude in degrees East
        :param lat: Latitude in degrees North
        :param grid: Set to false if lon,lat pairs represent individual points 
        :param engine: Spherical harmonic compute engine to use for the computation
        :return: A datarray for which the spherical harmonic coefficietn dimension is mapped to set of points
        The following scenarios can be handled:
        
        1: lon, lat are Xarray coordinate variables sharing the same dimension
        mension. Map to a list of points (SH dimension is mapped to a single dimension)
        2: lon, lat are Xarray coordinate variables with different dimensions: Map tot a grid spanned by lon,lat
        3 lon, lat are list-like objects with the same length: Map to a grid unless (grid=False)
        4. lon,lat are list-like objects of different lengths: Map to a grid
        """
        #dispatch to compute engine
        eng=self._eng(engine)()
        return eng.analysis(self._obj,lon,lat,grid)
    

@xr.register_dataset_accessor("sh")
class SHDsAccessor(ShXrBase):
    snrei=SnreiFactory
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    
    @staticmethod
    def zeros(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C'):
        """0-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,squeeze,name,auxcoords,order=order).to_dataset()
    
    @staticmethod
    def ones(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C'):
        """1-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,squeeze,name,auxcoords,order=order).to_dataset()
