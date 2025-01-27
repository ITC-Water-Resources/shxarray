# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

from shapely.geometry import Point
import xarray as xr
import numpy as np
import geopandas as gpd
import pandas as pd
from shxarray.kernels.axial import Disk,Unit,ParabolicCap

def point2sh(pointgeom,nmax:int=100,auxcoord=None,axialtype="unit",psi=None) ->xr.DataArray:
    """
        Convert a GeoSeries of points to axially symmetric loads expressed as SH coefficients
    Parameters
    ----------
    pointgeom : GeoSeries or array_like
       Iterable with points which describes the center of the loads
    auxcoord: named Pandas.Series or dict(dimname=coordvalues)
        Auxiliary coordinate to map to the dimension of pointgeom. The default will construct a coordinate with an sequential numerical index and index "id" 
    axialtype : str
        Type of the load to be constructed. One of 'unit','disk','paraboliccap'. Defaults to a unit load
    nmax : int
        maximum degree and order to resolve
        
    psi : int
        spherical width in degrees of the disk or parabolic cap. Ignored for unit loads. defaults to 1 degree

    Returns
    -------
    xr.DataArray
      A DataArray holding the spherical harmonic coefficients of the prescribed loads 
    
    See Also
    --------
    shxarray.geom.polygons.polygon2sh
    
    """
    
    if type(pointgeom) != gpd.GeoSeries:
        #convert to GeoSeries
        pointgeom=gpd.GeoSeries(pointgeom)

    if axialtype== "unit":
        axiso=Unit(nmax)
    elif axialtype == 'disk':
        axiso=Disk(nmax,psi)
    elif axialtype == "paraboliccap":
        axiso=ParabolicCap(nmax,psi)
    else:
        raise RuntimeError(f"Axial type {axialtype} is unknown")

    #Position the isotropic load on the chosen point locations
    daout=axiso.position(pointgeom.x,pointgeom.y)
    
    if auxcoord is None:
        dimk="id"
        coords=np.arange(len(pointgeom))
    else:
        if type(auxcoord) == pd.Series:
            dimk=auxcoord.name
            coords=auxcoord.values
        else:
            # should be a dictionary like object
            if len(auxcoord) != 1:
                raise RuntimeError("Only one input coordinate is accepted")
            dimk=next(iter(auxcoord))
            coords=auxcoord[dimk]
    return daout.assign_coords({dimk:("nlonlat",coords)})

     

