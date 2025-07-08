# Gravity functionals
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#





from shxarray.kernels.isokernelbase import IsoKernelBase,genCompIsoKernelClass
from shxarray.core.logging import shxlogger
from shxarray.earth.constants import a_earth,rho_water,rho_earth,rho_sea
from shxarray.earth.snrei import SnreiFactory
import xarray as xr
import numpy as np

gtypes=["stokes","tws","geoid","uplift","horzdef"]

class Stokes2TWS(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of Stokes coefficients [-] to equivalent water height [m]"""
    name="stokes2tws"
    transform=("stokes","tws")
    def __init__(self,knLove=None,nmax=None,k0=None,**kwargs):
        if knLove is None:
            #retrieve the default
            knLove=SnreiFactory.load(nmax=nmax,deg0scale=k0).kn
        shxlogger.warning("Thin shell approximation is used (surface load assumed)")
        super().__init__(dsiso=(2*knLove.n+1)*rho_earth*a_earth/(3*rho_water*(1+knLove)))
        
class Stokes2Geoid(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of disturbing potential to geoid height in meter, using Brun's formula"""
    name="stokes2geoid"
    transform=("stokes","geoid")
    def __init__(self,nmax,**kwargs):
        super().__init__(dsiso=a_earth*xr.DataArray(np.ones([nmax+1]),dims=['n'],coords=dict(n=np.arange(nmax+1))))

class TWS2Geoid(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of a surface load (in m) to geoid height in meter"""
    name="load2geoid"
    transform=("tws","geoid")
    def __init__(self,knLove=None,nmax=None,deg0scale=None,k0=None,**kwargs):
        if knLove is None:
            #retrieve the default

            knLove=SnreiFactory.load(nmax=nmax,deg0scale=k0).kn

        shxlogger.warning("Thin shell approximation is used (surface load assumed)")
        dsiso=(3*rho_water/rho_earth)*((1+knLove)/(2*knLove.n+1))
        if deg0scale is not None:
            dsiso.loc[0]=deg0scale
        super().__init__(dsiso=dsiso)

class TWS2Uplift(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of surface load (in m) to elastic uplift in meter"""
    name="tws2uplift"
    transform=("tws","uplift")
    def __init__(self,hnLove=None,nmax=None,h0=None,**kwargs):
        if hnLove is None:
            #retrieve the default
            hnLove=SnreiFactory.load(nmax=nmax,deg0scale=h0).hn
        shxlogger.warning("Thin shell approximation is used (surface load assumed)")
        super().__init__(dsiso=(3*rho_water/rho_earth)*(hnLove/(2*hnLove.n+1)))

class TWS2Horzdef(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of surface load (in m) to elastic horizontal deformation component
    Note to compute the horizontal deformation, the spatial derivatives in longitude and latitude need to be computed
    """
    name="tws2horzdef"
    transform=("tws","horzdef")
    def __init__(self,lnLove=None,nmax=None,l0=None,**kwargs):
        if lnLove is None:
            #retrieve the default
            lnLove=SnreiFactory.load(nmax=nmax,deg0scale=l0).ln
        shxlogger.warning("Thin shell approximation is used (surface load assumed)")
        super().__init__(dsiso=(3*rho_water/rho_earth)*(lnLove/(2*lnLove.n+1)))


#create the inverse classes of the above
TWS2Stokes=Stokes2TWS.invcls()
Geoid2Stokes=Stokes2Geoid.invcls()
Uplift2TWS=TWS2Uplift.invcls()
Geoid2TWS=TWS2Geoid.invcls()
Horzdef2TWS=TWS2Horzdef.invcls()

#generate some composite classes
Stokes2Uplift=genCompIsoKernelClass(Stokes2TWS,TWS2Uplift)
Uplift2Stokes=Stokes2Uplift.invcls()
Geoid2Uplift=genCompIsoKernelClass(Geoid2TWS,TWS2Uplift)
Geoid2Horzdef=genCompIsoKernelClass(Geoid2TWS,TWS2Horzdef)
Horzdef2Geoid=Geoid2Horzdef.invcls()
Stokes2Horzdef=genCompIsoKernelClass(Stokes2TWS,TWS2Horzdef)
Horzdef2Stokes=Stokes2Horzdef.invcls()

gravclasses=[Stokes2TWS,TWS2Geoid,TWS2Uplift,Stokes2Geoid,TWS2Stokes,Geoid2Stokes,Uplift2TWS,Geoid2TWS,Stokes2Uplift,Geoid2Uplift,TWS2Horzdef,Horzdef2TWS,Geoid2Horzdef,Uplift2Stokes,Stokes2Horzdef,Horzdef2Stokes,Horzdef2Geoid]

gravlookup={cls.transform:cls for cls in gravclasses}


def gravFunc(fromType,toType,**kwargs):
    """Computes a kernel to transform of one gravitational function into another"""
    transtype=(fromType,toType)
    if transtype in gravlookup:
        return gravlookup[transtype](**kwargs)

    raise RuntimeError(f"No gravity transfer function available for transform {fromType}->{toType}")

