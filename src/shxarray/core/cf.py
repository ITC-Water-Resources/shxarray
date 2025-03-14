# Contain some basic functionality for finding and using CF variables
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import numpy as np
from collections import namedtuple
from shxarray.core.admin import get_userinfo
from shxarray._version import version

CoordInfo=namedtuple("CoordInfo",['min','max','step','direction','var'])

def find_coord(coordvars, names):
    """
        Find a coordinate variable in a dictionary of coordinate variables
    """
    tmp=[ky for ky in coordvars.keys() if ky in names]
    if len(tmp) > 1 or len(tmp) == 0:
        raise KeyError(f"Cannot find an unambigious coordinate variable from {names}")
    coordvar=coordvars[tmp[0]]
    
    #also compute minmax and possibly uniform step size
    cmin=coordvar.min().item()
    cmax=coordvar.max().item()
    
    cincr=np.diff(coordvar.data)
    dmin=cincr.min()
    dmax=cincr.max()
    if dmin == dmax:
        #uniform step
        cstep=dmin
        if dmin < 0:
            direction='descending'
        else:
            direction='ascending'
    else:
        #nonuniform step
        cstep=None
        if dmin < 0 and dmax < 0:
            direction='descending'
        elif dmin > 0 and dmax > 0:
            direction='ascending'
        else:
            direction="random"

    return CoordInfo(min=cmin,max=cmax,step=cstep,direction=direction,var=coordvar)



def find_lon(coordvars):
    return find_coord(coordvars,['lon','Longitude','x','longitude'])

def find_lat(coordvars):
    return find_coord(coordvars,['lat','Latitude','y','latitude'])

def change_central_longitude(dain,central_longitude=0,resort=True):
    """
    Change the central longitude of a dataset to 180 or 360 degrees
    Parameters
    ----------
    dain : xarry.DataArray or xarray.Dataset
        input data to change
    central_longitude : int
        central longitude to change to. either 0 or 180 degrees (default is 0)
        

    Returns
    -------
        xarray.DataArray or xarray.Dataset
            Input data with changed central longitude (if needed) otherwise returns the original input

    """
    if central_longitude not in [180,0]:
        raise ValueError("Central longitude should be 180 or 0")

    loninfo = find_lon(dain.coords)
    #guess for correct central longitude
    if loninfo.min < 0 or loninfo.max <= 180:
        given_central_longitude=0
    else:
        given_central_longitude=180
    if given_central_longitude != central_longitude:
        if central_longitude == 180:
            dain.coords[loninfo.var.name]=(dain.coords[loninfo.var.name] + 360) % 360
        else:
            dain.coords[loninfo.var.name]=(dain.coords[loninfo.var.name] + 180) % 360 - 180
        if resort:
            dain=dain.sortby(dain.coords[loninfo.var.name])
    
    return dain


cflookup={
        "longitude":{'units':'degrees_east','standard_name':'longitude','long_name':'longitude'},
        "latitude":{'units':'degrees_north','standard_name':'latitude','long_name':'latitude'},
        "stokes":{'units':'-',"long_name":"Stokes Coefficients","gravtype":"stokes"},
        "stokes stdv":{'units':'-',"long_name":"Standard deviation of the Stokes Coefficients","gravtype":"stokes"},
        "tws":{'units':'m',"long_name":"Total water storage","gravtype":"tws"},
        
        }

def get_cfatts(standard_name):
    """Return CF attributes for certain coordinate types"""
    return cflookup[standard_name]

def get_cfglobal():
    """Return global attributes and possible user details for the CF convention"""
    cfattr={'Conventions':'CF-1.8',"source":f"shxarray-{version()} <https://github.com/ITC-Water-Resources/shxarray>"} 
    cfattr.update(get_userinfo())
    return cfattr



