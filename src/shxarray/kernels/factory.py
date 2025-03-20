# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

from shxarray.kernels.ddk import load_ddk
from shxarray.kernels.gauss import Gaussian

def getSHfilter(filtername,nmax,**kwargs):
    """
    Retrieve the (aniso) tropic kernel for known SH filters
    Parameters
    ----------
    filtername : str
        e.g. Gauss200 for a Gaussian filter with a halfwidth of 200 km, DDK5 for a DDK filter with sharpness 5
        
    nmax : int
        maximum degree of the filter to accomodate
        
    **kwargs : dict
        additional keyword arguments to pass to the filter
        

    Returns
    -------
    A kernel filter operator

    """
    if filtername.startswith('DDK'):
        #load a dedicated DDK filter
        if "transpose" in kwargs:
            trans=kwargs["transpose"]
        else:
            trans=False
        if "truncate" in kwargs:
            truncate=kwargs["truncate"]
        else:
            truncate=True
        kernel=load_ddk(filtername,trans,nmax,truncate)
    elif filtername.startswith('Gauss'):
        if "halfwidth" in kwargs:
            radius=kwargs["halfwidth"]
        else:
            try:
                radius=int(filtername[5:])
            except:
                raise RuntimeError("Cannot parse the Gaussian halfwidth in km.\n Specify either 'Gaussxxx'or add halfwidth=xxx to the sh.filter call")
        kernel=Gaussian(nmax,radius*1e3)
    else:
        raise RuntimeError(f"SH Filter {filtername} not recognized")
    
    return kernel
