# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#





import xarray as xr

def multiply(dash1:xr.DataArray,dash2:xr.DataArray,engine="shtns")->xr.DataArray:
    """
    Multiply two spherical harmonics DataArrays together (equivalent to multiplying in the spatial domain)
    Currently this function uses a spherical harmonic synthesis followed, a multiplication in the spatial domain and a spherical harmonic analysis
    A future implementation may use the realGaunt coefficients to directly erform the multiplication in the spectral domain
    Parameters
    ----------
        
    dash1 : xr.DataArray
        
    dash2 : xr.DataArray
        
    engine : str
        Engine to use for the synthesis and analysis step
    

    Returns
    -------
    xr.DataArray
        

    """
    dagrd1=dash1.sh.synthesis(engine=engine)
    dagrd2=dash2.sh.synthesis(engine=engine)
    dagrd=dagrd1*dagrd2
    dashout=dagrd.sh.analysis(nmax=dash1.sh.nmax,engine=engine)
    return dashout


