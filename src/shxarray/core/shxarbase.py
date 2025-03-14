# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import xarray as xr
import numpy as np

from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.logging import shxlogger

loaded_engines={}

# version which is backward compatible:

try:
    from importlib.metadata import entry_points
except:
    breakpoint()
    from importlib_metadata import entry_points

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
        
        if "n" in self._obj.coords:
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
       
        if "n" in self._obj.coords:
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
    
    @gravtype.setter
    def gravtype(self,gravtypeval):
        """Sets the gravitational type of the content"""
        if gravtypeval not in ["stokes","tws"]:
            shxlogger.warning(f"Unknown gravitation type {gravtypeval}")
        self._obj.attrs["gravtype"]=gravtypeval
   
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
    def _initWithScalar(nmax,nmin=0,scalar=0,name="cnm",auxcoords={},order='C',nshdims=1):
        """Initialize an spherical harmonic DataArray based on nmax and nmin"""
        if nshdims <= 0 or nshdims > 2:
            raise RuntimeError("Spherical harmonic coefficient dimension must be 1 or 2")
        coords={SHindexBase.name:SHindexBase.nm_mi(nmax,nmin)}
        if nshdims == 2:
            #add an additional SH coordinate
            coords[SHindexBase.name_t]=SHindexBase.mi_toggle(coords[SHindexBase.name]) 
        dims=[]
        shp=[]
        
        #possibly use auxiliary coordinates
        for dim,coord in auxcoords.items():
            dims.append(dim)
            shp.append(len(coord))
            coords[dim]=coord
        
        # add shi dimension and shape last (so it varies quikest in memory
        if nshdims == 2:
            dims.append(SHindexBase.name_t)
            shp.append(len(coords[SHindexBase.name_t]))
        dims.append(SHindexBase.name)
        shp.append(len(coords[SHindexBase.name]))


        if scalar == 0:
            return xr.DataArray(data=np.zeros(shp,order=order),dims=dims,name=name,coords=coords)
        elif scalar == 1:
            return xr.DataArray(data=np.ones(shp,order=order),dims=dims,name=name,coords=coords)
        else:
            return xr.DataArray(data=np.full(shp,scalar,order=order),dims=dims,name=name,coords=coords)

    


    
    def drop_nmindex(self,suf=''):
        if suf:
            indname=f"{SHindexBase.name}{suf}"
        else:
            indname=SHindexBase.name
        return self._obj.reset_index(indname)
    
    def build_nmindex(self,suf=''):
        if suf:
            indname=f"{SHindexBase.name}{suf}"
            nname=f"n{suf}"
            mname=f"m{suf}"

        else:
            indname=SHindexBase.name
            nname='n'
            mname='m'

        if indname in self._obj.indexes:
            #already build, so don't bother
            return self._obj
        
        #either build from separate coordinate variables (n,m,t)
        if nname in self._obj.coords and mname in self._obj.coords:
            shimi=SHindexBase.mi_fromtuples([(n,m) for n,m in zip(self._obj[nname].values,self._obj[mname].values)],suf)
            if hasattr(xr,'Coordinates'):
                #only in newer xarray versions..
                shimi=xr.Coordinates.from_pandas_multiindex(shimi, indname)
            else:
                shimi={indname:(indname,shimi)}

            return self._obj.drop_vars([nname,mname]).assign_coords(shimi)

        elif indname in self._obj.coords:
            #rebuild multiindex from an array of "left-over" tuples
            shimi=SHindexBase.mi_fromtuples(self._obj[indname].values,suf)
            if hasattr(xr,'Coordinates'):
                #only in newer xarray versions..
                shimi=xr.Coordinates.from_pandas_multiindex(shimi, indname)
            else:
                shimi={indname:(indname,shimi)}
            return self._obj.drop_vars([indname]).assign_coords(shimi)
    
    def set_nmindex(self,shimi,suf=''):
        """Sets the spherical harmonic coordinate indesd from a given pandas multiindex""" 
        if suf:
            indname=f"{SHindexBase.name}{suf}"
            nname=f"n{suf}"
            mname=f"m{suf}"

        else:
            indname=SHindexBase.name
            nname='n'
            mname='m'
        
        if shimi.names[0] != nname or shimi.names[1] != mname:
            raise RuntimeError(f"Level names of index should be named:{nname},{mname}")
        if hasattr(xr,'Coordinates'):
            #only in newer xarray versions..
            shimi=xr.Coordinates.from_pandas_multiindex(shimi, indname)
        else:
            shimi={indname:(indname,shimi)}

        return self._obj.assign_coords(shimi)

        

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
            
            if engine == "exp":
                eng=__import__("shxarray.exp")
                loaded_engines[engine]=eng.exp
            else:
                grp="shxarray.computebackends"
                eps=entry_points(group=grp)

                if engine not in eps.names:
                    raise RuntimeError(f"Compute engine {engine} not found, is it installed?")
                Eng=eps[engine].load()
                loaded_engines[engine]=Eng()
        


        return loaded_engines[engine]
    
    @staticmethod
    def lonlat_grid(nmax,engine="shlib",**kwargs):
        """
        Return a dataset with the required lon/lat coordinates needed for a given maximum degree

        Parameters
        ----------
        nmax : int
            Maximum degree of the spherical harmonic expansion
            
        engine : str, optional
            The compute engine to use. The default is "shlib".
        **kwargs : dict
            Additional keyword arguments to be passed to the engine. For example, the grid type can be set with the gtype keyword

        Returns
        -------
        xarray.Dataset
            A xarray.Dataset with the required lon/lat span as cooridnates

        """
        eng=ShXrBase._eng(engine)
        return eng.lonlat_grid(nmax,**kwargs)

