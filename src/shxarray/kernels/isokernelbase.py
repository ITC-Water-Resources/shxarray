# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.logging import logger
from scipy.sparse import diags
from shxarray.shlib import Ynm

class IsoKernelBase:
    """
    Provides functionality to work with isotropic spherical harmonic kernels
    """
    attr={"shtype":"shiso","kernelstate":"collapsed"}
    name="shkernel"
    def __init__(self):
        self._dsiso=None
    
    @property
    def nmax(self):
        return self._dsiso.n.max().item()
    
    @property
    def nmin(self):
        return self._dsiso.n.min().item()
   
    def expanddiag(self,shindex):
        nmin=shindex.n.min().item()
        nmax=shindex.n.max().item()
        nminsup= self._dsiso.n.min().item() 
        nmaxsup= self._dsiso.n.max().item()
        if nmin < nminsup or  nmax > nmaxsup:
            raise RuntimeError("Requested kernel operation is only supported for degrees {nminsup} < = n <= {nmaxsup}")

        if self._dsiso.n.diff(dim="n").max().item() > 1:
            logger.info("Some degrees are missing in the kernel, interpolating")
            coeff=self._dsiso.interp(n=shindex.n)
        else:
            coeff=self._dsiso.sel(n=shindex.n)
    
        return xr.DataArray(coeff.data,coords=dict(shi=shindex))


    def jacobian(self,shindex):
        # create a sparse diagnonal array
        return diags(self.expanddiag(shindex).values)

    def __call__(self,dain:xr.DataArray):
        #create the jacobian matrix based on the input maximum and minimum degrees
        if "shi" not in dain.indexes:
            raise RuntimeError("Spherical harmonic index not found in input, cannot apply kernel operator to object")
        #expand kernel to the same degrees as the input
        daexpand=self.expanddiag(dain.shi)
        daout=dain*daexpand
        return daout

    def position(self,lon,lat):
        """
        Position this kernel on a specific location of the sphere
        :param lon: Longitude(s) in degrees of position (list like)
        :param lat: Latitude(s) in degrees of position (list-like)
        :return: A xarray.DataArray with the kernel located on the specified locations
        """
        ynm=Ynm(self.nmax)
        ynmdata=ynm(lon,lat)
        # scale elements by 1/(2n+1)
        normv=[2*n+1 for n in ynmdata.n.data]
        ynmdata=ynmdata/normv

        return self.__call__(ynmdata)

    def Greensfunc(self,theta):
        """
        Map an isotropic kernel to a 1D Greens function (as a function of distance from the center)
        :param theta: discretized isotropic distance in degrees
        :return: A xarray.DataArray with the 1-D Greens function in the spatial domain
        """
        pass
