# Axial symmetric expansions of unit and disk loads
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


from shxarray.kernels.isokernelbase import IsoKernelBase
import xarray as xr
import numpy as np
from shxarray.shlib import Pn
from math import cos, radians
class Unit(IsoKernelBase):
    name="shunit"
    def __init__(self,nmax):
        """
        Create an isotropic kernel representing a unit load
        :param nmax: maximum degree to resolve
        :return: A Isotropic Kernel object representing a unit load
        """
        super().__init__()
        scales=[2*n+1 for n in range(nmax+1)]
        self._dsiso=xr.DataArray(scales,coords={"n":("n",range(nmax+1))},dims=["n"])

class Disk(IsoKernelBase):
    name="shdisk"
    def __init__(self,nmax,psi):
        """
        Create an isotropic kernel representing a disk load
        :param nmax: maximum degree to resolve
        :param psi: disk size in angular degrees
        :return: A Isotropic Kernel object representing a disk load
        """
        super().__init__()
        self.psi=psi
        cospsi=np.cos(radians(psi))
        legendre=Pn(nmax+1)(cospsi)
        scales=np.zeros([nmax+1])
        scales[0]=(1-cospsi)/2
        for n in range(1,nmax+1):
            scales[n]=(legendre[n-1]-legendre[n+1])/2

        self._dsiso=xr.DataArray(scales,coords={"n":("n",range(nmax+1))},dims=["n"])

class ParabolicCap(IsoKernelBase):
    name="shcap"
    def __init__(self,nmax,psi):
        """
        Create an isotropic kernel representing a Parabolic Cap
        :param nmax: maximum degree to resolve
        :param psi: Cap size in angular degrees
        :return: A Isotropic Kernel object representing a parabolic cap
        """
        super().__init__()
        self.psi=psi
        psi=radians(psi)
        cospsi=cos(psi)
        scales=np.zeros([nmax+1])
        scales[0]=(1-cospsi)/3
        for n in range(1,nmax+1):
            scales[n]=((cos((n-1)*psi)-cos(n*psi))/(n-0.5)-(cos((n+1)*psi)-cos((n+2)*psi))/(n+1.5))/((1.0-cos(psi))*4.0)

        self._dsiso=xr.DataArray(scales,coords={"n":("n",range(nmax+1))},dims=["n"])


