# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.logging import logger
from shxarray.shlib import Ynm
# from dask.array.core import einsum_lookup
import sparse




class AnisoKernel:
    """
    Provides functionality to work with anisotropic spherical harmonic kernels
    """
    attr={"shtype":"shaniso"}
    name="shanisokernel"
    def __init__(self,dsobj):
        self._dskernel=dsobj
    
    @property
    def nmax(self):
        return self._dskernel.sh.nmax
    
    @property
    def nmin(self):
        return self._dskernel.sh.nmin
   

    def __call__(self,dain:xr.DataArray):
        if "shi" not in dain.indexes:
            raise RuntimeError("al harmonic index not found in input, cannot apply kernel operator to object")
        # daout=xr.dot(dain,self._dskernel.mat,dims=["shi"]) 
        daout=xr.dot(self._dskernel.mat,dain,dims=["shi"]) 
        #rename shi and convert to dense array
        daout=daout.sh.toggle_shi()
        daout=xr.DataArray(daout.data.todense(),coords=daout.coords)
        return daout

    def position(self,lon,lat):
        """
        Position this kernel on a specific location of the sphere
        :param lon: Longitude(s) in degrees of position (list like)
        :param lat: Latitude(s) in degrees of position (list-like)
        :return: A xarray.DataArray with the kernel located on the specified locations
        """
        ynm=Ynm(self.nmax)
        ynmdata=ynm(lon,lat)
        # scale elements by 1/(2n+1)
        normv=[2*n+1 for n in ynmdata.n.data]
        ynmdata=ynmdata/normv

        return self.__call__(ynmdata)

