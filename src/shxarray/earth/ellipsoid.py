
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

import xarray as xr
import shxarray

def get_GRS80_Stokes():
    """
    Return the GRS80 ellipsoid as Stokes coefficients
    """
    a = 6378137.0
    f = 1.0 / 298.257222101
    nmax=8
    dsgrs80=xr.DataArray.sh.zeros(nmax)
    dsgrs80.attrs['title']='GRS80'
    dsgrs80.attrs['a_earth']=a
    dsgrs80.attrs['f_earth']=f
    #fill out derived values
    dsgrs80.loc[dict(n=0,m=0)]=1
    dsgrs80.loc[dict(n=2,m=0)]=-0.48416685489612e-03 
    dsgrs80.loc[dict(n=4,m=0)]=0.79030407333333e-06 
    dsgrs80.loc[dict(n=6,m=0)]=-0.16872510013651e-08 
    dsgrs80.loc[dict(n=8,m=0)]=-0.34609833692685e-11 

    return dsgrs80
    
