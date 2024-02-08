# Smoothing kernels
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


from shxarray.kernels.isokernelbase import IsoKernelBase
import numpy as np
from shxarray.earth.constants import a_earth
from math import cos,log,exp
import xarray as xr
class Gaussian(IsoKernelBase):
    name="gauss"
    def __init__(self,nmax,halfwidth):
        """
        Create an isotropic kernel representing a unit load. According to Wahr et al 1998 with a cutoff criteria to avoid numerical instability
        :param nmax: maximum degree to resolve
        :param halfwidth: Halfwidth of the Gaussian kernel in m (at the Earth's surface)
        :return: A Isotropic Kernel object representing a unit load
        """
        super().__init__()
        self.halfwidth=halfwidth
        arg=log(2.0)/(1-cos(halfwidth/a_earth))
        exparg=exp(-2*arg)
        
        wn=np.zeros([nmax+1])
        wn[0]=1
        wn[1]=(1+exparg)/(1-exparg) -1/arg
        for n in range(1,nmax):
            wn[n+1]= wn[n-1] -(2*n+1)/arg*wn[n]
            if wn[n+1] < 1e-8:
               break
        self._dsiso=xr.DataArray(wn,coords={"n":("n",range(nmax+1))},dims=["n"])
