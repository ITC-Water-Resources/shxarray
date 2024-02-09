# Contain some basic functionality for finding and using CF variables
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import numpy as np
from collections import namedtuple

CoordInfo=namedtuple("CoordInfo",['min','max','step','var'])

def find_coord(coordvars, names):
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
        else:
            #nonuniform step
            cstep=None

        return CoordInfo(min=cmin,max=cmax,step=cstep,var=coordvar)



def find_lon(coordvars):
    return find_coord(coordvars,['lon','Longitude','x','longitude'])

def find_lat(coordvars):
    return find_coord(coordvars,['lat','Latitude','y','latitude'])

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

