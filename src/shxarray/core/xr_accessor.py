# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import xarray as xr
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.shxarbase import ShXrBase
from shxarray.core.logging import shxlogger
from shxarray.kernels import getSHfilter
from shxarray.io.shascii import to_shascii
from shxarray.geom.polygons import polygon2sh
from shxarray.geom.points import point2sh

import numpy as np
from shxarray.kernels.gravfunctionals import gravFunc,gtypes
from shxarray.signal.basinav import Basinav

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
    
    def synthesis(self,lon=None, lat=None,engine="shlib",gtype=None,**kwargs):
        """

        Apply spherical harmonic synthesis on a set of longitude, latitude points
        gen
        :param lon: Longitude in degrees East
        :param lat: Latitude in degrees North
        :param gtype: Set to false if lon,lat pairs represent individual points 
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
        if lon is None and lat is None:
            #automatically find a suitable grifd based on the input maximum degree
            dslonlat=eng.lonlat_grid(nmax=self._obj.sh.nmax,gtype=gtype)
        elif lon is not None and lat is not None:
            #generate a grid consistent with the backend
            if len(lon) == len(lat) and gtype is None:
                gtype='point'
            dslonlat=eng.lonlat_grid(lon=lon,lat=lat,gtype=gtype)
        else:
            raise RuntimeError("Both or neither of lon and lat should be specified")

        return eng.synthesis(self._obj,dslonlat,**kwargs)
    
    def synthesis_like(self,dslonlat=None,engine="shlib",**kwargs):
        #dispatch to compute engine
        eng=self._eng(engine)
        return eng.synthesis(self._obj,dslonlat,**kwargs)

    def analysis(self,nmax=100,engine="shlib",**kwargs):
        """
        Apply spherical harmonic analysis from the given ints
        :param nmax : Spherical harmonic truncation degree of output
        :return: A datarray with spherical harmonic coefficients derived from the input dataArray

        Depending on the method applied different scenarios can be handled
        input is given on an equidistant longitude, latitude grid (but may be different step size in lon and lat direction)
        """
        #dispatch to compute engine
        eng=self._eng(engine)
    
        return eng.analysis(self._obj,nmax=nmax,**kwargs)

    def filter(self,filtername,**kwargs):
        """
        Apply well known filters to Spherical harmonic data
        :param filtername: currently 'DDKX' or 'Gauss'
        :param **kwargs:
            transpose (default=False): apply he transpose of the filter (only makes sense for anisotropic filters
            halfwidth (int): specify the halfwidth in km's for the Gaussian filter, alternatively specify in the format Gauss300
            truncate (bool) (default=True): Truncate low degree coefficients which fall outside the filter. Set to False to keep unfiltered input
        """
        kernel=getSHfilter(filtername,self.nmax,**kwargs)
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
        return self.gravfunctional("tws",**kwargs)
    
    def geoid(self,**kwargs):
        return self.gravfunctional("geoid",**kwargs)
    
    def uplift(self,**kwargs):
        return self.gravfunctional("uplift",**kwargs)
    
    def horzdef(self,**kwargs):
        return self.gravfunctional("horzdef",**kwargs)
    
    def stokes(self,**kwargs):
        return self.gravfunctional("stokes",**kwargs)

    def gravfunctional(self,outgravtype,ingravtype=None,**kwargs):
        """
        Compute a gravity functional from the spherical harmonic coefficients
        :param outgravtype: The gravity functional type to compute (see shxarray.kernels.gravfunctionals.gtypes)
        :param ingravtype: The input gravity functional type to use (default is the one set in the DataArray)
        :return: A DataArray with the computed gravity functional
        """
        if ingravtype is None:
            #try retrieving the input gravitational type
            try:
                ingravtype=self.gravtype
            except TypeError:
                raise TypeError(f"Input gravtype cannot be retrieved, user setter xarray.DataArray.sh.gravtype=xxx or add ingravtype=xxx as an argument, where  xxx is one of ({','.join(gtypes)})")
        if ingravtype == outgravtype:
            shxlogger.warning(f"Input and output gravity functional type are the same ({ingravtype}), returning input")
            return self._obj

        if self._obj.sh.nmin == 0:
            if "tws" in outgravtype or "tws" in ingravtype or "uplift" in outgravtype or "uplift" in ingravtype or "horzdef" in outgravtype or "horzdef" in ingravtype:
                

                #set the k0 value to 0.0 if not specified
                #this is needed for TWS and uplift computations
                if "k0" not in kwargs:
                    kwargs["k0"]=0.0
                    shxlogger.warning('k0 Load Love number is not set explicitly, but n=0 is requested, setting k0=0.0 for Total water storage computation (treat degree 0 output with caution)')
                if "h0" not in kwargs:
                    #set the h0 value to 1.0 if not specified
                    kwargs["h0"]=1.0
                    shxlogger.warning('h0 Load Love number is not set explicitly, but n=0 is requested, setting h0=1.0 for uplift computation (treat degree 0 output with caution)')

                if "l0" not in kwargs:
                    #set the h0 value to 1.0 if not specified
                    kwargs["l0"]=1.0
                    shxlogger.warning('l0 Load Love number is not set explicitly, but n=0 is requested, setting l0=1.0 for horizontal deformation computation, but treat degree 0 output with caution')

        kernel=gravFunc(ingravtype,outgravtype,nmax=self.nmax,**kwargs)
        return self.convolve(kernel)


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
                
    def dvplot(self,ax=None,mean=False,**kwargs):
        """
            Plot the degree variance of spherical harmonic coefficients

        Parameters
        ----------
        ax : 
            Matplotlib axis object to plot on. If None, a new axis will be created
            
        **kwargs:
            Additional arguments passed to the matplotlib function
            
        Returns:
        --------
            ax: Matplotlib axis object
        """

        dadv=self.degvar(mean)
        if ax is None:
            lplt=dadv.plot(**kwargs)
        else:
            lplt=dadv.plot(ax=ax,**kwargs)
        if ax is None:
            ax=lplt.axes
        fig=lplt.figure
        return ax

    def geoplot(self,ax=None,engine="shlib",add_colorbar=True,**kwargs):
        """
            Plot the spherical harmonic coefficients on a global map using cartopy
        Parameters
        ----------
        ax :
            Cartopy axis object to plot on. If None, a new axis will be created
        engine : str, default: 'shlib'
            Compute engine to use for the synthesis of the spherical harmonic coefficients. Other options could be 'shtns' (when installed).
        add_colorbar : bool, default: True
            Whether to add a colorbar to the plot
        **kwargs:
            Additional arguments passed to the pcolormesh plot function
        Returns
        -------
        ax : cartopy axis object
            The axis on which the spherical harmonic coefficients are plotted
        """

        import cartopy.crs as ccrs
        # import cartopy.feature as cfeature
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.axes(projection=ccrs.PlateCarree())



        dagrd=self._obj.sh.synthesis(engine=engine)

        if add_colorbar:
            cbar_kwargs = dict(orientation="horizontal",label=f"{self._obj.name} [{self._obj.attrs.get('units','?')}]")
        else:
            cbar_kwargs = None

        dagrd.plot.pcolormesh(
            transform=ccrs.PlateCarree(), add_colorbar=add_colorbar,ax=ax,cbar_kwargs=cbar_kwargs,**kwargs)
        ax.set_global()
        if 'long_name' in self._obj.attrs:
            ax.set_title(self._obj.attrs['long_name'])
        ax.coastlines(color='#282828', linewidth=0.5)
        
        return ax

    @staticmethod    
    def from_geoseries(gseries,nmax:int,auxcoord=None,engine="shlib",**kwargs):
        """
            Convert a GeoSeries (from geopandas) to spherical harmonic coefficients
        Parameters
        ----------
        gseries : geopandas.GeoSeries 
            A GeoSeries Instance of points or polygons
            
        nmax : int
            maximum spherical harmonic degree and order to resolve
            
        auxcoord : named Pandas.Series or dict(dimname=coordvalues)
            Auxiliary coordinate to map to the dimension of gseries. The default will construct a coordinate with an sequential numerical index and index "id" 
        engine: str, default: 'shlib'
            Compute engine to use for the computation. Other options could be 'shtns' (when installed). This option has no effect when the input is a GeoSeries of points
        **kwargs: dict
           Optional arguments which will be passed to either polygon2sh, or point2sh   

        Returns
        -------
        xr.DataArray
        A DataArray holding the spherical harmonic coefficients up to maximum degree specified

        See Also
        --------
        shxarray.geom.points.point2sh
        shxarray.geom.polygons.polygon2sh

        """
        
        gtypes=gseries.geom_type.unique()
        if len(gtypes) > 1:
            raise RuntimeError("from_gseries does not currently accept mixed geometry types")
        if gtypes[0] == "Polygon" or gtypes[0] == "MultiPolygon":
            return polygon2sh(gseries,nmax=nmax,auxcoord=auxcoord,engine=engine,kwargs=kwargs)

        elif gtypes[0] =="Point":
            return point2sh(gseries,nmax=nmax,auxcoord=auxcoord,**kwargs)
        else:
            raise RuntimeError(f"geometry type {gtypes[0]}, is not supported")

    def basinav(self,dabasins,filtername=None,leakage_corr=None,**kwargs):
        """
        Compute the basin averages of the spherical harmonic coefficients using basin coefficients and possibly a filter and leakage correction
        Parameters
        ----------
        dabasins : xr.DataArray 
            spherical harmonic basin coefficients (basins are uniform masks with 1 inside and 0 outside the region)
        filtername : str or None
            Name of the filter to apply to the coefficients (see sh.filter)
            
        leakage_corr : str or None
            Name of the leakage correction to apply to the coefficients.
            

        """
        basinavOp=Basinav(dabasins,filtername=filtername,leakage_corr=leakage_corr)
        return basinavOp(self._obj,**kwargs)
   
    def multiply(self,daother,engine="exp",truncate=False,**kwargs):
        #dispatch to compute engine
        eng=self._eng(engine)
        if truncate:
            nmax=self._obj.sh.nmax
            daout=eng.multiply(self._obj,daother,**kwargs)
            return daout.sh.truncate(nmax=nmax)
        else:
            return eng.multiply(self._obj,daother,**kwargs)

@xr.register_dataset_accessor("sh")
class SHDsAccessor(ShXrBase):
    def __init__(self, xarray_obj):
        super().__init__(xarray_obj)
    
    def synthesis(self,lon=None, lat=None,engine="shlib",**kwargs):
        """Calls the spherical harmonic synthesis operation on all DataArrays which have a 'nm' index"""

        #gather relevant das
        das={ky:da for ky,da in self._obj.data_vars.items() if SHindexBase.name in da.coords}
        dsout=None
        for name,da in das.items():
            daout=da.sh.synthesis(lon=lon, lat=lat,engine=engine,**kwargs)
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
