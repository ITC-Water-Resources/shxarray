# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import xarray as xr
import numpy as np

from .sh_indexing import SHindexBase


loaded_engines={}

# version which is backward compatible:
from importlib_metadata import entry_points
# For newer versions this may eventually need to be replaced with 
# from importlib.metadata import entry_points

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
        if SHindexBase.name in self._obj.indexes:
           return self._obj.nm.n.min().item()

        if "n" in self._obj.indexes:
            return self._obj.n.min().item()
        
        raise RuntimeError("Cannot return nmin, spherical harmonic index is not initialized") 
    
    @property
    def nmax(self):
        """
        maximum degree of the spherical harmonic dataset
        :return: maximum degree
        :rtype: int
        """
        if SHindexBase.name in self._obj.indexes:
            return self._obj.nm.n.max().item()

        if "n" in self._obj.indexes:
            return self._obj.n.max().item()
        
        raise RuntimeError("Cannot return nmax, spherical harmonic index is not initialized") 
    
    @property
    def gravtype(self):
        """Returns the registered gravitational type of the content"""
        try:
            return self._obj.attrs["gravtype"]
        except KeyError:
            # defaults to Stokes coefficients 
            return "stokes"
   
    def truncate(self,nmax=None,nmin=None,dims=[SHindexBase.name]):
        """
        Truncate the maximum and/or minimum degree of the spherical harmonic coordinate and corresponding variables
        :param nmax: (int) maximum spherical harmonic degree to keep
        :param nmin: (int) minimum spherical harmonic degree to keep
        :param dims: list of SH dimensions to truncate over default ['shi']. Alternatively specify ['shi_'] or both
        :return: A new truncated data array
        :rtype: xarray.DataArray
        """
        indx=None
        da=None
        if SHindexBase.name in dims:
            if nmax is not None:
                indx=(self._obj.nm.n <= nmax)
            if nmin is not None:
                if indx is not None:
                    indx=indx*(self._obj.nm.n >= nmin)
                else:
                    indx=(self._obj.nm.n >= nmin)
            da=self._obj.isel(nm=indx)
        
        if "nm_" in dims:
            if nmax is not None:
                indx=(self._obj.nm_.n_ <= nmax)
            if nmin is not None:
                if indx is not None:
                    indx=indx*(self._obj.nm_.n_ >= nmin)
                else:
                    indx=(self._obj.nm_.n_ >= nmin)
            if da is None:
                da=self._obj.isel(nm_=indx)
            else:
                da=da.isel(nm_=indx)
         
        if da is not None:
            return da
        else:
            raise RuntimeError("No spherical harmonic index ('nm' or 'n') was found in the xarray object")
    
    @staticmethod
    def _initWithScalar(nmax,nmin=0,scalar=0,name="cnm",auxcoords={},order='C'):
        """Initialize an spherical harmonic DataArray based on nmax and nmin"""
        
        coords={SHindexBase.name:SHindexBase.nm_mi(nmax,nmin)}
        dims=[]
        shp=[]
        
        
        #possibly use auxiliary coordinates
        for dim,coord in auxcoords.items():
            dims.append(dim)
            shp.append(len(coord))
            coords[dim]=coord
        
        # add shi dimension and shape last (so it varies quikest in memory
        dims.append(SHindexBase.name)
        shp.append(len(coords[SHindexBase.name]))
        
        if scalar == 0:
            return xr.DataArray(data=np.zeros(shp,order=order),dims=dims,name=name,coords=coords)
        elif scalar == 1:
            return xr.DataArray(data=np.ones(shp,order=order),dims=dims,name=name,coords=coords)
        else:
            return xr.DataArray(data=np.full(shp,scalar,order=order),dims=dims,name=name,coords=coords)

    
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


    
    def drop_nmindex(self):
        ds=self._obj.reset_index(SHindexBase.name)
        return ds.assign_coords(t=([SHindexBase.name],[t for t in ds.t.values]))
    
    def build_nmindex(self):
        if SHindexBase.name in self._obj.indexes:
            #already build, so don't bother
            return self._obj
        #either build from separate coordinate variables (n,m,t)
        if "n" in self._obj.coords and "m" in self._obj.coords and "t" in self._obj.coords:
            shimi=SHindexBase.mi_fromtuples([(n,m) for n,m in zip(self._obj.n.values,self._obj.m.values)])
            return self._obj.drop_vars(["n","m"]).assign_coords(nm=shimi)
        elif SHindexBase.name in self._obj.coords:
            #rebuild multiindex from an array of "left-over" tuples
            shimi=SHindexBase.mi_fromtuples(self._obj.nm.values)
            return self._obj.drop_vars([SHindexBase.name]).assign_coords(nm=shimi)
    
    def toggle_nm(self):
        """Toggle naming of nm, nm_ multindices and their levels"""
        renamedict={}
        if SHindexBase.name in self._obj.dims:
            renamedict[SHindexBase.name]=SHindexBase.name_t
            renamedict["n"]="n_"
            renamedict["m"]="m_"
        
        if SHindexBase.name_t in self._obj.dims:
            renamedict[SHindexBase.name_t]=SHindexBase.name
            renamedict["n_"]="n"
            renamedict["m_"]="m"
        return self._obj.rename(renamedict)

    @staticmethod
    def _eng(engine="shlib"):
        

        if engine not in loaded_engines:
            #load engine if not done already

            grp="shxarray.computebackends"
            eps=entry_points(group=grp)

            if engine not in eps.names:
                raise RuntimeError(f"compute engine {engine} not found")
            Eng=eps[engine].load()
            loaded_engines[engine]=Eng()
        


        return loaded_engines[engine]
        

