# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import xarray as xr
import numpy as np

from .sh_indexing import SHindexBase,trig
from importlib.metadata import entry_points


class ShXrBase:
    def __init__(self, xarray_obj):
        """Constructor takes an xarray object as conforming with the xarray accessor routines
        :param xarray_obj: Xarray.DataArray  representing the parent
        """
        self._obj = xarray_obj
    
    @property
    def nmin(self):
        """
        Minimum degree of the spherical harmonic data
        :return: minimum degree
        :rtype: int
        """
        if self._obj.shi is None:
            raise RuntimeError("Cannot return nmin, spherical harmonic index is not initialized") 
        return self._obj.shi.n.min().item()
    
    @property
    def nmax(self):
        """
        maximum degree of the spherical harmonic dataset
        :return: maximum degree
        :rtype: int
        """
        if self._obj.shi is None:
            raise RuntimeError("Cannot return nmax, spherical harmonic index is not initialized") 
        return self._obj.shi.n.max().item()
   
    def truncate(self,nmax=None,nmin=None):
        """
        Truncate the maximum and/or minimum degree of the spherical harmonic coordinate and corresponding variables
        :param nmax: (int) maximum spherical harmonic degree to keep
        :param nmin: (int) minimum spherical harmonic degree to keep
        :return: A new truncated data array
        :rtype: xarray.DataArray
        """
        indx=None
        if nmax is not None:
            indx=(self._obj.shi.n <= nmax)
        if nmin is not None:
            if indx is not None:
                indx=indx*(self._obj.shi.n >= nmin)
            else:
                indx=(self._obj.shi.n >= nmin)
        
        return self._obj.isel(shi=indx) 


    @staticmethod
    def _initWithScalar(nmax,nmin=0,scalar=0,squeeze=True,name="cnm",auxcoords={}):
        """Initialize an spherical harmonic DataArray based on nmax and nmin"""
        
        coords={"shi":SHindexBase.nmt_mi(nmax,nmin,squeeze=squeeze)}
        dims=["shi"]
        shp=[len(coords['shi'])]

        #possibly append coordinates
        for dim,coord in auxcoords.items():
            dims.append(dim)
            shp.append(len(coord))
            coords[dim]=coord


        if scalar == 0:
            return xr.DataArray(data=np.zeros(shp),dims=dims,name=name,coords=coords)
        elif scalar == 1:
            return xr.DataArray(data=np.ones(shp),dims=dims,name=name,coords=coords)
        else:
            return xr.DataArray(data=np.full(shp,scalar),dims=dims,name=name,coords=coords)

    
    # @staticmethod
    # def from_cnm(cnm):
        # """Create a xarray from a cnm array from shtools"""
        # nmax=cnm.shape[1]-1
        # #create a multiindex
        # shgmi=SHAccessor.nmt_mi(nmax,squeeze=True)
        # #indexing vectors so the sh coeffients match the index
        # i_n=shgmi.get_level_values(level='n').astype(int)
        # i_m=shgmi.get_level_values(level='m').astype(int)
        # i_t=shgmi.get_level_values(level='t').astype(int)
        # coords={"shg":shgmi}
        # dims=["shg"]
        # return xr.DataArray(data=cnm[i_t,i_n,i_m],dims=dims,name="cnm",coords=coords)


    
    def drop_shindex(self):
        ds=self._obj.reset_index("shi")
        return ds.assign_coords(t=(["shi"],[t for t in ds.t.values]))
    
    def build_shindex(self):
        if "shi" in self._obj.indexes:
            #already build, so don't bother
            return self._obj
        #either build from separate coordinate variables (n,m,t)
        if "n" in self._obj.coords and "m" in self._obj.coords and "t" in self._obj.coords:
            shimi=SHindexBase.mi_fromtuples([(n,m,trig(t)) for n,m,t in zip(self._obj.n.values,self._obj.m.values,self._obj.t.values)])
            return self._obj.drop_vars(["n","m","t"]).assign_coords(shi=shimi)
        elif "shi" in self._obj.coords:
            #rebuild multiindex from an array of "left-over" tuples
            shimi=SHindexBase.mi_fromtuples(self._obj.shi.values)
            return self._obj.drop_vars(["shi"]).assign_coords(shi=shimi)

    @staticmethod
    def _eng(engine="shlib"):
        computebackends=entry_points(group="shxarray.computebackends",name=engine)
        if len(computebackends) == 0:
            raise RuntimeError(f"compute engine {engine} not found")
        return computebackends[0].load()

    
