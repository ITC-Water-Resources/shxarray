# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.shxarbase import ShXrBase
from shxarray.kernels.ddk import load_ddk
from shxarray.kernels.gauss import Gaussian
from shxarray.io.shascii import to_shascii
from shxarray.geom.polygons import polygon2sh
from shxarray.geom.points import point2sh

import numpy as np
from shxarray.kernels.gravfunctionals import gravFunc

@xr.register_dataarray_accessor("sh")
class SHDaAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    

    @staticmethod
    def zeros(nmax,nmin=0,name="cnm",auxcoords={},order='C',nshdims=1):
        """0-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,name,auxcoords,order=order,nshdims=nshdims)
    
    @staticmethod
    def ones(nmax,nmin=0,name="cnm",auxcoords={},order='C',nshdims=1):
        """1-Initialize an spherical harmonic DataArray based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,name,auxcoords,order=order,nshdims=nshdims)
    
    @staticmethod
    def wigner3j(j2,j3,m2,m3, engine="shlib"):
        #dispatch to compute engine
        eng=ShXrBase._eng(engine)
        return eng.wigner3j(j2,j3,m2,m3)
    
    @staticmethod
    def gauntReal(n2,n3,m2,m3, engine="shlib"):
        #dispatch to compute engine
        eng=ShXrBase._eng(engine)
        return eng.gauntReal(n2,n3,m2,m3)
    
    @staticmethod
    def gaunt(n2,n3,m2,m3, engine="shlib"):
        #dispatch to compute engine
        eng=ShXrBase._eng(engine)
        return eng.gaunt(n2,n3,m2,m3)
    
    def synthesis(self,lon=None, lat=None,grid=True,engine="shlib"):
        """
        Apply spherical harmonic synthesis on a set of longitude, latitude points
        :param lon: Longitude in degrees East
        :param lat: Latitude in degrees North
        :param grid: Set to false if lon,lat pairs represent individual points 
        :param engine: Spherical harmonic compute engine to use for the computation
        :return: A datarray for which the spherical harmonic coefficient dimension is mapped to set of points
        The following scenarios can be handled:
        
        1: lon, lat are Xarray coordinate variables sharing the same dimension
        mension. Map to a list of points (SH dimension is mapped to a single dimension)
        2: lon, lat are Xarray coordinate variables with different dimensions: Map tot a grid spanned by lon,lat
        3 lon, lat are list-like objects with the same length: Map to a grid unless (grid=False)
        4. lon,lat are list-like objects of different lengths: Map to a grid
        """
        #dispatch to compute engine
        eng=self._eng(engine)
         
        if lon is None:
            lon=np.arange(-180.0,181.0,1.0)
        if lat is None:
            lat=np.arange(-90.0,91.0,1.0)
        return eng.synthesis(self._obj,lon,lat,grid)
    
    def analysis(self,nmax=100,method='integrate',engine="shlib"):
        """
        Apply spherical harmonic analysis from the given ints
        :param nmax : Spherical harmonic truncation degree of output
        :param method: Method to use for the analysis
        :return: A datarray with spherical harmonic coefficients derived from the input dataArray

        Depending on the method applied different scenarios can be handled
        method == 'integrate'
        input is given on an equidistant longitude, latitude grid (but may be different step size in lon and lat direction)
        """
        #dispatch to compute engine
        eng=self._eng(engine)
    
        return eng.analysis(self._obj,nmax,method)

    def filter(self,filtername,**kwargs):
        """
        Apply well known filters to Spherical harmonic data
        :param filtername: currently 'DDKX' or 'Gauss'
        :param **kwargs:
            transpose (default=False): apply he transpose of the filter (only makes sense for anisotropic filters
            halfwidth (int): specify the halfwidth in km's for the Gaussian filter, alternatively specify in the format Gauss300
            truncate (bool) (default=True): Truncate low degree coefficients which fall outside the filter. Set to False to keep unfiltered input
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
            kernel=load_ddk(filtername,trans,self.nmax,truncate)
        elif filtername.startswith('Gauss'):
            if "halfwidth" in kwargs:
                radius=kwargs["halfwidth"]
            else:
                try:
                    radius=int(filtername[5:])
                except:
                    raise RuntimeError("Cannot parse the Gaussian halfwidth in km.\n Specify either 'Gaussxxx'or add halfwidth=xxx to the sh.filter call")
            kernel=Gaussian(self.nmax,radius*1e3)
        else:
            raise RuntimeError(f"SH Filter {filtername} not recognized")

        return self.convolve(kernel)

    def convolve(self,kernel):
        """
        Execute a convolution of the dataarray with a given kernel/filter:
        :param kernel: Kernel object (see  e.g shxarray.kernels)
        """
        return kernel(self._obj)

    def degvar(self,mean=False):
        """
        Compute the degree variances of spherical harmonic data
        :param mean: Take the average power per degree instead of the sum (divide by 2n+1)
        :return: A Dataarray with the degree variances
        """

        if mean:
            dv=np.square(self._obj).sh.drop_nmindex().set_xindex("n").groupby("n").mean()
            dv.attrs=self._obj.attrs 
            dv.attrs["long_name"]="Average degree variance"
        else:
            dv=np.square(self._obj).sh.drop_nmindex().set_xindex("n").groupby("n").sum()
            dv.attrs=self._obj.attrs 
            dv.attrs["long_name"]="Degree variance"
        
        dv.n.attrs["long_name"]="degrees (n)"
        dv.name='cn'
        return dv 
    
    def tws(self,**kwargs):
        kernel=gravFunc(self.gravtype,"tws",nmax=self.nmax,**kwargs)
        return self.convolve(kernel)
    
    # def geoid(self,**kwargs):
        # pass #return self.gravfunctional("geoid",**kwargs)

    def to_ascii(self,out_obj=None):
        """Return a string representing the ascii file content of the spherical harmonic coefficients""" 
        return to_shascii(self._obj,out_obj)
    
    def p2s(self,engine="shlib"):
        """Returns the product to sum matrix of the spherical harmonic coefficients in the dataarray
        The output matrix will be symmetric with sides spanning up to nmax/2 of the input"""
        eng=self._eng(engine)
        return eng.p2s(self._obj.sh.build_nmindex())
    
    def triplot(self,ax=None,**kwargs):
        """
            Plot the spherical harmonic coefficients as a triangular plot

        Parameters
        ----------
        ax : 
            Matplotlib axis object to plot on. If None, a new axis will be created
            
        **kwargs:
            Additional arguments passed to the pcolormesh plot function
            
        Returns:
        --------
            ax: Matplotlib axis object
        """

        qdmesh=self._obj.unstack('nm').plot(ax=ax,add_colorbar=False,**kwargs)
        if ax is None:
            ax=qdmesh.axes
        ax.set_aspect('equal')
        fig=qdmesh.figure
        fig.colorbar(qdmesh,orientation='horizontal')
        return ax
                
    @staticmethod    
    def from_geoseries(gseries,nmax:int,auxcoord=None,**kwargs):
        """
            Convert a GeoSeries (from geopandas) to spherical harmonics
        Parameters
        ----------
        gseries : geopandas.GeoSeries 
            A GeoSeries Instance 
            
        nmax : int
            maximum spherical harmonic degree and order to resolve
            
        auxcoord : named Pandas.Series or dict(dimname=coordvalues)
            Auxiliary coordinate to map to the dimension of gseries. The default will construct a coordinate with an sequential numerical index and index "id" 
        
        **kwargs: dict
           Optional arguments which will be passed to either polygon2sh, or point2sh   

        Returns
        -------
        xr.DataArray
        A DataArray holding the spherical harmonic coefficients

        """
        
        gtypes=gseries.geom_type.unique()
        if len(gtypes) > 1:
            raise RuntimeError("from_gseries does not currently accept mixed geometry types")
        if gtypes[0] == "Polygon" or gtypes[0] == "MultiPolygon":
            return polygon2sh(gseries,nmax=nmax,auxcoord=auxcoord,kwargs=kwargs)
        elif gtypes[0] =="Point":
            return point2sh(gseries,nmax=nmax,auxcoord=auxcoord,**kwargs)
        else:
            raise RuntimeError(f"geometry type {gtypes[0]}, is not supported")

        import pdb;pdb.set_trace() 
        pass
    

@xr.register_dataset_accessor("sh")
class SHDsAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    
    def synthesis(self,lon=None, lat=None,grid=True,engine="shlib"):
        """Calls the spherical harmonic synthesis operation on all DataArrays which have a 'nm' index"""
        

        #gather relevant das
        das={ky:da for ky,da in self._obj.data_vars.items() if SHindexBase.name in da.coords}
        dsout=None
        for name,da in das.items():
            daout=da.sh.synthesis(lon=lon, lat=lat,grid=grid,engine=engine)
            if dsout is None:
                dsout=daout.to_dataset(name=name)
            else:
                dsout[name]=daout
        return dsout
    
    @staticmethod
    def zeros(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C',nshdims=1):
        """0-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,0,squeeze,name,auxcoords,order=order,nshdims=nshdims).to_dataset()
    
    @staticmethod
    def ones(nmax,nmin=0,squeeze=True,name="cnm",auxcoords={},order='C',nshdims=1):
        """1-Initialize an spherical harmonic Dataset based on nmax and nmin"""
        return ShXrBase._initWithScalar(nmax,nmin,1,squeeze,name,auxcoords,order=order,nshdims=nshdims).to_dataset()
