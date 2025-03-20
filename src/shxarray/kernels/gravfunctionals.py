# Gravity functionals
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#





from shxarray.kernels.isokernelbase import IsoKernelBase

from shxarray.earth.constants import a_earth,rho_water,rho_earth,rho_sea
from shxarray.earth.snrei import SnreiFactory
import xarray as xr
import numpy as np

class Stokes2TWS(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of Stokes coefficients [-] to equivalent water height [m]"""
    name="stokes2tws"
    transform=("stokes","tws")
    def __init__(self,knLove=None,nmax=None):
        super().__init__()
        if knLove is None:
            #retrieve the default
            knLove=SnreiFactory.load(nmax=nmax).kn

        self._dsiso=(2*knLove.n+1)*rho_earth*a_earth/(3*rho_water*(1+knLove))

class Stokes2Geoid(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of disturbing potential to geoid height in meter, using Brun's formula"""
    name="stoked2geoid"
    transform=("stokes","geoid")
    def __init__(self,nmax):
        super().__init__()
        self._dsiso=a_earth*xr.DataArray(np.ones([nmax+1]),dims=['n'],coords=dict(n=np.arange(nmax+1)))

class Load2Geoid(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of a surface load (in m) to geoid height in meter"""
    name="load2geoid"
    transform=("load","geoid")
    def __init__(self,knLove=None,nmax=None,deg0scale=None):
        super().__init__()
        if knLove is None:
            #retrieve the default
            knLove=SnreiFactory.load(nmax=nmax,deg0scale=deg0scale).kn

        self._dsiso=(3*rho_water/rho_earth)*((1+knLove)/(2*knLove.n+1))
        if deg0scale is not None:
            self._dsiso.loc[0]=deg0scale

class Load2Uplift(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of surface load (in m) to elastic uplift in meter"""
    name="load2uplift"
    transform=("load","uplift")
    def __init__(self,hnLove=None,nmax=None,deg0scale=None):
        super().__init__()
        if hnLove is None:
            #retrieve the default
            hnLove=SnreiFactory.load(nmax=nmax,deg0scale=deg0scale).hn
        self._dsiso=(3*rho_water/rho_earth)*(hnLove/(2*hnLove.n+1))
        
gravclasses=[Stokes2TWS,Load2Geoid,Load2Uplift,Stokes2Geoid]

gravlookup={cls.transform:cls for cls in gravclasses}

def gravFunc(fromType,toType,**kwargs):
    """Computes a kernel to transform of one gravitational function into another"""
    transtype=(fromType,toType)
    if transtype in gravlookup:
        return gravlookup[transtype](**kwargs)
    elif transtype[::-1] in gravlookup:
        #inverse version is found
        return gravlookup[transtype[::-1]](**kwargs).inv()

    raise RuntimeError(f"No gravity transfer function available for transform {fromType}->{toType}")

