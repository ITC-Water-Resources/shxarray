# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.core.logging import shxlogger
from scipy.sparse import diags
from shxarray.shlib import Ynm
from shxarray.core.cf import get_cfatts
from shxarray.core.sh_indexing import SHindexBase

class IsoKernelBase:
    """
    Provides functionality to work with isotropic spherical harmonic kernels
    """
    attr={"shtype":"shiso","kernelstate":"collapsed"}
    name="shkernel"
    transform=None
    invert=False #inverts the kernel when initializing
    def __init__(self,dsiso=None,name=None,transform=None):
        self._dsiso=dsiso
        if self.invert:
            if self._dsiso is None:
                raise ValueError("Cannot invert kernel without a dsiso defined")
            self._dsiso=1/self._dsiso

        if name is not None:
            self.name=name
        if transform is not None:
            self.transform=transform
    
    @classmethod
    def invcls(cls):
        """
        Returns the inverse class of the current kernel class
        :return: An IsoKernelBase class representing the inverse kernel
        """
        if cls.transform is None:
            trans=None
        else:
            trans= cls.transform[::-1]
        clsnew=type(f"{cls.__name__}_inv",(cls,),dict(invert=True,transform=trans,name=f"{cls.name}_inv"))
        return clsnew

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
            raise RuntimeError(f"Requested kernel operation is only supported for degrees {nminsup} < = n <= {nmaxsup}")

        if self._dsiso.n.diff(dim="n").max().item() > 1:
            shxlogger.info("Some degrees are missing in the kernel, interpolating")
            coeff=self._dsiso.interp(n=shindex.n)
        else:
            coeff=self._dsiso.sel(n=shindex.n)
         
        return xr.DataArray(coeff.data,coords=dict(nm=shindex))


    def jacobian(self,shindex):
        # create a sparse diagnonal array
        return diags(self.expanddiag(shindex).values)

    def __call__(self,dain:xr.DataArray):
        #create the jacobian matrix based on the input maximum and minimum degrees
        if SHindexBase.name not in dain.indexes:
            raise RuntimeError("Spherical harmonic index not found in input, cannot apply kernel operator to object")
        #expand kernel to the same degrees as the input
        daexpand=self.expanddiag(dain.nm)
        daout=dain*daexpand
        if self.transform is not None:
            name=self.transform[1]

        else:
            name=self.name
        try:
            #try to update the dataarray attributes
            daout.attrs.update(get_cfatts(name))
        except:
            pass
        return daout.rename(name)
    
    def inv(self):
        """Returns the inverse of an already instantiated isotropic Kernel"""
        invkernel=self.invcls()(self._dsiso)
        return invkernel


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

    def plot(self,ax=None,**kwargs):
        """
        Plot the isotropic kernel as a function of degree

        Parameters
        ----------
        ax : matplotlib axis object, optional
            The axis to plot the kernel on. If None, a new axis is created.
            
        **kwargs: additional keyword arguments for plotting
            

        Returns
        -------
        ax : matplotlib axis object
            
        """
        
        lplt=self._dsiso.plot(ax=ax,**kwargs)
        if ax is None:
            ax=lplt[0].axes

        ax.set_title(f"{self.name} Kernel")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Kernel Coefficient")
        return ax




class IsoKernelCompositeBase(IsoKernelBase):
    """
    Combines two isotropic kernels into a joint one
    """
    
    #note clsleft and clsright need to be defined in the derived class
    clsleft=None
    clsright=None
    def __init__(self, **kwargs):
        

        if self.clsleft is None or self.clsright is None:
            raise ValueError("A derived composite kernel classes needs to have clsleft and clsright defined")

        # initialize the left and right kernels
        # Note: both kernels need to accept the same arguments (but not all need to be used)
        lft=self.clsleft(**kwargs)
        rght=self.clsright(**kwargs)
        
        dsisonew= lft._dsiso * rght._dsiso
        name= f"{rght.transform[1]} -> {lft.transform[0]}"
        super().__init__(dsiso=dsisonew,name=name)


def genCompIsoKernelClass(clsleft,clsright):
    """
    Create a composite kernel  class (not instance!) from two isotropic kernels
    :param clsleft: The left IsoKernelBase class to use
    :param clsright: The right IsoKernelBase class to use
    :return: An IsoKernelCompositeBase object representing the composite kernel
    """
    

    if clsleft.transform is None or clsright.transform is None:
        raise ValueError("Cannot create composite kernel without a transform defined for both kernels")

    

    if clsleft.transform[1] != clsright.transform[0]:
        raise ValueError(f"Compositekernel join dimensions are inconsistent")
    
    transform=(clsleft.transform[0],clsright.transform[1])

    return type(f"Comp_{clsleft.__name__}_{clsright.__name__}", (IsoKernelCompositeBase,), dict(clsleft=clsleft,clsright=clsright,transform=transform))

