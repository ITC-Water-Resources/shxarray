# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#





import xarray as xr
from shxarray.shlib import getGauntReal
from shxarray.core.logging import shxlogger
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import numpy as np


def multiply_spat(dash1:xr.DataArray,dash2:xr.DataArray,engine="shtns")->xr.DataArray:
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

def multiply(dash1:xr.DataArray,dash2:xr.DataArray):
    
    norm=np.sqrt(4*np.pi)
    nm1=dash1.nm
    nm2=dash2.nm
    # breakpoint() 
    #allocate space for the output
    nmaxout=dash1.sh.nmax+dash2.sh.nmax
    auxcoords=dash1.sh.auxcoords()
    dashout=xr.DataArray.sh.zeros(nmax=nmaxout,auxcoords=auxcoords)
    with logging_redirect_tqdm():
        for n1,m1 in tqdm(nm1.data):
            shxlogger.info(f"Multiplying n1,m1: {n1},{m1}")
            for n2,m2 in nm2.data:
                gnt = getGauntReal(n1,n2,m1,m2)
                if 28 in gnt.n and 28 in gnt.m: 
                    breakpoint()
                dashout.loc[{"nm":gnt.nm}]+= norm*gnt*dash1.loc[{"nm":(n1,m1)}]*dash2.loc[{"nm":(n2,m2)}]
    return dashout

