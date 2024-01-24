# Gravity functionals
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#





from shxarray.kernels.isokernelbase import IsoKernelBase

from shxarray.earth.constants import a_earth,rho_water,rho_earth,rho_sea

class Stokes2Eqh(IsoKernelBase):
    """Provides an isotropic kernel representing the transformation of Stokes coefficients [-] to equivalent water height [m]"""
    name="stokes2eqh" 
    def __init__(self,knLove):
        super().__init__()
        self._dsiso=(2*knLove.n+1)*rho_earth*a_earth/(3*rho_water*(1+knLove))

