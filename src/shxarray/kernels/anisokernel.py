# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.core.logging import shxlogger
from shxarray.shlib import Ynm
from shxarray.core.sh_indexing import SHindexBase
import sparse

from packaging import version





class AnisoKernel:
    """
    Provides functionality to work with anisotropic spherical harmonic kernels
    """
    attr={"shtype":"shaniso"}
    def __init__(self,dsobj,name="aniso",truncate=True):
        self._dskernel=dsobj
        self.name=name 
        self.truncate=truncate
        self.useDask=version.parse(xr.__version__) < version.parse('2023.11.0')
        if self.useDask:
            from dask.array.core import einsum_lookup
            #Register the einsum functions which are needed to do the sparse dot functions (for earlier versions of xarray)
            einsum_lookup.register(sparse.COO,AnisoKernel.daskeinsumReplace)
            #convert to xarray with dask structure
            self._dskernel=self._dskernel.chunk()

    @property
    def nmax(self):
        return self._dskernel.sh.nmax
    
    @property
    def nmin(self):
        return self._dskernel.sh.nmin
   

    def __call__(self,dain:xr.DataArray):
        if SHindexBase.name not in dain.indexes:
            
            raise RuntimeError("SH harmonic index not found in input, cannot apply kernel operator to object")

        if self.nmax < dain.sh.nmax:
            raise RuntimeError("Input data has higher degree than kernel, cannot apply kernel operator to object")

        daout=xr.dot(self._dskernel.mat,dain,dims=[SHindexBase.name]) 
        #rename nm and convert to dense array
        daout=daout.sh.toggle_nm()
        if self.useDask:
            daout=xr.DataArray(daout.compute().data,coords=daout.coords,name=self.name)
        else:
            if hasattr(daout.data,'todense'):
                #still needs expanding
                daout=xr.DataArray(daout.data.todense(),coords=daout.coords,name=self.name)
            else:
                #just rename
                daout.name=self.name

        if not self.truncate and self.nmin > dain.sh.nmin:
            #also add the unfiltered lower degree coefficients back to the results
            daout=xr.concat([dain.sh.truncate(self.nmin-1,dain.sh.nmin),daout],dim=SHindexBase.name)

        if daout.sh.nmax == dain.sh.nmax and daout.sh.nmin == dain.sh.nmin:
            #resort to original order when output nmax and nmin agree
            daout=daout.sel({SHindexBase.name:dain.nm})

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
    
    @staticmethod
    def daskeinsumReplace(subscripts, *operands, out=None, dtype=None, order='K', casting='safe', optimize=False):
        """Mimics the interface of https://numpy.org/doc/stable/reference/generated/numpy.einsum.html, but uses the sparse.COO dot function"""
        if subscripts == "ab,cb->ac":
            return operands[0].dot(operands[1].T)
        elif subscripts == "ab,ca->bc":
            return operands[0].T.dot(operands[1].T)
        elif subscripts == "ab,bc->ac":
            return operands[0].dot(operands[1])
        elif subscripts == "ab,b->a":
            return operands[0].dot(operands[1])
        elif subscripts == "ab,a->b":
            return operands[0].T.dot(operands[1])
        elif subscripts == "ab,ac->bc":
            return operands[0].T.dot(operands[1])
          
        else:
            raise NotImplementedError(f"Don't know (yet) how to handle this einsum: {subscripts} with sparse.dot operations")

