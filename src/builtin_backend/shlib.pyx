# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++
# cython: profile=False
"""
**shlib** is shxarray's default binary Cython backend. 
Some of the heavy lifting such as synthesis and analysis operations, is done using this the functions of this shared library.
"""

include "legendre.pyx" 
include "wigner3j.pyx"
include "gaunt.pyx"
include "ynm.pyx"
include "synthesis.pyx"
include "analysis.pyx"
include "sinex.pyx"

from shxarray.core.shcomputebase import SHComputeBackendBase
from shxarray.core.cf import get_cfatts,get_cfglobal
import numpy as np

class SHComputeBackend(SHComputeBackendBase):
    _credits="Used backend: shlib cython extension from shxarray (https://shxarray.wobbly.earth/latest/references/shxarray.html#module-shxarray.shlib)"

    
    def synthesis(self,dain,dslonlat,**kwargs):
        
        syn=Synthesis(dslonlat)
        synres=syn(dain)
        synres.name=dain.name
        synres.lon.attrs=get_cfatts("longitude")
        synres.lat.attrs=get_cfatts("latitude")
        synres.attrs.update(get_cfglobal())
        synres.attrs['history']="Synthesis performed using shlib backend"
        synres.attrs['comments']=self._credits
        return synres

    def analysis(self,dain,nmax,**kwargs):
        ana=Analysis(nmax)
        anares=ana(dain)
        anares.attrs.update(get_cfglobal())
        anares.attrs['history']="Analysis performed using shlib backend"
        anares.attrs['comments']=self._credits
        return anares
        
    
    def gaunt(self,n2,n3,m2,m3):
        return getGaunt(n2,n3,m2,m3)

    def gauntReal(self,n2,n3,m2,m3):
        return getGauntReal(n2,n3,m2,m3)

    def wigner3j(self,j2,j3,m2,m3):
        return getWigner3j(j2,j3,m2,m3)

    def p2s(self,daobj):
        return getp2s(daobj)
    
    def lonlat_grid(self,nmax=None,lon=None,lat=None,gtype="regular"):
        """
        Create a lon-lat grid compatible with the SHlib backend
        Parameters
        ----------
        nmax : int, optional
            Maximum expected degree of the spherical harmonic expansion.
            
        gtype : str, optional
            Type of grid to create. Only 'regular','regular_lon0' and "point" are supported. The regular_lon0 option is a regular grid with the first longitude point at 0 degrees.
        Returns
        -------
            An empty xarray.Dataset with CF-complying lon, lat coordinates variables

        """
        if gtype is None:
            gtype="regular"


        if gtype not in ["regular","regular_lon0","point"]:
            raise ValueError("Only 'regular','regular_lon0' and 'point' type is supported")
        
        if nmax is not None:
            if lon is not None or lat is not None:
                raise ValueError("If nmax is specified, lon and lat should not be specified")
            #autogenerate a grid based on nmax
            idres=1
            while idres > 360/nmax/4:
                idres=idres/2
            if gtype == "regular":
                lon=np.arange(-180+idres/2,180,idres)
                lat=np.arange(-90+idres/2,90,idres)
                
            elif gtype == "regular_lon0":
                lon=np.arange(0,360,idres)
                lat=np.arange(-90+idres/2,90,idres)
            else:
                raise ValueError("Only 'regular' and 'regular_lon0' grid type are currently supported")

            dslonlat=xr.Dataset(coords=dict(lon=lon,lat=lat))
        elif lon is not None and lat is not None:
            if gtype == "point":
                if len(lon) != len(lat):
                    raise ValueError("For point type 'grid', lon and lat should have the same length")
                dslonlat=xr.Dataset(coords=dict(lon=("nlonlat",lon.data),lat=("nlonlat",lat.data)))
            else:
                dslonlat=xr.Dataset(coords=dict(lon=lon,lat=lat))
            dslonlat=xr.Dataset(coords=dict(lon=lon,lat=lat))
        else:
            raise ValueError("Either nmax or lon and lat should be specified")
        
        dslonlat.lon.attrs.update(get_cfatts("longitude"))
        dslonlat.lat.attrs.update(get_cfatts("latitude"))

        dslonlat.lon.attrs["shxarray_gtype"]=gtype
        dslonlat.lat.attrs["shxarray_gtype"]=gtype
        return dslonlat


    def multiply(self,dash1:xr.DataArray,dash2:xr.DataArray,method="spectral")->xr.DataArray:
        """
        Multiply two spherical harmonics DataArrays together (equivalent to multiplying in the spatial domain)
        This function uses the computation of Real Gaunt coefficients
        Parameters
        ----------
            
        dash1 : xr.DataArray
            
        dash2 : xr.DataArray
            


        Returns
        -------
        xr.DataArray
            

        """
        from time import time
        if method == "spatial":
            # Note: approximate method only
            shxlogger.info("Multiplying in the spatial domain")
            nmax=dash1.sh.nmax+dash2.sh.nmax
            dslonlat=self.lonlat_grid(nmax=nmax)
            t0=time()
            dagrd1=self.synthesis(dash1,dslonlat=dslonlat)
            dagrd2=self.synthesis(dash2,dslonlat=dslonlat)
            t1=time()
            dagrd=dagrd1*dagrd2
            dashout=self.analysis(dagrd,nmax=nmax)
            t2= time()
            shxlogger.info(f"Multiplication in the spatial domain took {t1-t0:.2f} seconds for synthesis and {t2-t1:.2f} seconds for analysis")
        elif method == "spectral":
            shxlogger.info("Multiplying in the spectral domain")
            dashout=multiply_sh(dash1,dash2)
        else:
            raise ValueError("Method for shlib should be either 'spatial' or 'spectral'")

        return dashout

