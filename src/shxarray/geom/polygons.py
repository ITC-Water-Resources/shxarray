# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

from shapely.geometry import Point
import xarray as xr
import shxarray
import numpy as np
import geopandas as gpd
import pandas as pd

from shxarray.core.logging import logger

def polygon2sh(polygeom,nmax:int=100,auxcoord=None) ->xr.DataArray:
    """
    Convert a mask defined by a polygon to spherical harmonic coefficients.
    This routine currently uses a simple integration approach

    
    Parameters
    ----------
    polygeom : GeoSeries or array_like
       Iterable with polygons which describes the mask (1 inside/0 outside)
       If the input is a GeoSeries, the crs of the provided geoseries will be used to do the polygon test. E.g. use a Antarctic Polar Stereographic epsg 3031, to do test for polgyons which contain the South Pole
        
    nmax : int
        Maximum degree and order of the output
    auxcoord: named Pandas.Series or dict(dimname=coordvalues)
        Auxiliary coordinate to map to the dimension of polygeom. The default will construct a coordinate with an sequential numerical index and index "id" 
    Returns
    -------
    xr.DataArray
        A DataArray holding the spherical harmonic coefficients
        

    
    """
    
    #create a dense enough grid encompassing all polgyons to use for spherical harmonic synthesis
    # heuristic way to figure out resolution based on nmax
    idres=1
    while idres > 360/nmax/4:
        idres=idres/2
    
    if type(polygeom) != gpd.GeoSeries:
        polygeom=gpd.GeoSeries(polygeom)

    #create a grid
    lon=np.arange(-180+idres/2,180,idres)
    lat=np.arange(-90+idres/2,90,idres)
    dims=["lon","lat"]
    coords={"lon":lon,"lat":lat}
    
    if auxcoord is None:
        coords["id"]=np.arange(len(polygeom))
        dims.append("id")
    else:
        if type(auxcoord) == pd.Series:
            dimk=auxcoord.name
            coords[dimk]=auxcoord.values
        else:
            # should be a dictionary like object
            if len(auxcoord) != 1:
                raise RuntimeError("Only one input coordinate is accepted")
            dimk=next(iter(auxcoord))
            coords[dimk]=auxcoord[dimk]
            
        dims.append(dimk)


    dtmp=xr.DataArray(np.zeros([len(lon),len(lat),len(polygeom)]),coords=coords,dims=dims).stack(lonlat=("lon","lat"))

    #create a geoDataframe of points from the grid
    ggrd=gpd.GeoDataFrame(geometry=[Point(lon,lat) for lon,lat in dtmp.lonlat.values],crs=4326)
    
    if polygeom.crs != ggrd.crs:
        #possibly convert the lon/lat grid in the desired projection before doing the polygon test
        ggrd=ggrd.to_crs(polygeom.crs)
    

    #query using a spatial index and set values to 1
    logger.info("Gridding polygons")
    for i,poly in enumerate(polygeom): 
        idx=ggrd.sindex.query(poly,predicate="contains")
        dtmp[i,idx]=1.0
    
    dtmp=dtmp.unstack("lonlat")
    # import ipdb;ipdb.set_trace() 
    logger.info("Applying SH analysis")
    dsout=dtmp.sh.analysis(nmax) 

    return dsout






