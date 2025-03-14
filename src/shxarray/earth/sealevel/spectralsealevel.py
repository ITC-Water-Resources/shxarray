# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#


from shxarray.earth.sealevel.sealevel import SeaLevelSolver
import xarray as xr
from math import floor
import os
from shxarray.core.admin import defaultcache
from shxarray.core.logging import shxlogger
from shxarray.kernels.gravfunctionals import Load2Geoid,Load2Uplift
from shxarray.earth.rotation import qs_rotfeedback_slow


class SpectralSeaLevelSolver(SeaLevelSolver):

    def __init__(self,oceansh:xr.DataArray=None, nmax=None,dssnrei=None,p2scache=None,rotfeedback=False):
        super().__init__(rotfeedback)
        
        if oceansh is None:
            if nmax is None or nmax > 150:
                raise RuntimeError("nmax (<=150) must be provided when no oceansh datatarray is provided")
            #use default ocean function (download)
            oceanshfile=os.path.join(defaultcache("ocean"),"ne_10m_oceansh_n300.nc")
            if not os.path.exists(oceanshfile):
                import requests
                url="https://github.com/strawpants/geoshapes/raw/refs/heads/master/ocean/ne_10m_oceansh_n300.nc"
                r = requests.get(url)
                shxlogger.info(f"Downloading ocean SH coefficients {oceanshfile}")
                with open(oceanshfile,'wb') as fid:
                    fid.write(r.content)
            else:
                shxlogger.info(f"{oceanshfile}, already downloaded")
            
            oceansh=xr.open_dataset(oceanshfile).sh.truncate(nmax=2*nmax).oceansh
        else:
            if nmax is not None:
                if nmax > oceansh.sh.nmax/2:
                    raise RuntimeError(f"Requested maximum degree {nmax} not supported by file (needs to be at least 2*nmax)")
                oceansh=oceansh.sh.truncate(nmax=nmax*2)
    
        if nmax is None:
            # Note: for full spectral consistency, the maximum degree of the ocean function is half that of the input ocean function
            self.nmax=floor(oceansh.sh.nmax/2)
        else:
            self.nmax=nmax
       
        #note: ocean coefficients need to be supported up to 2*nmax
        
        #setup SNREI Earth loading function
        if dssnrei is None:
            #default uses PREM Earth Model
            self.geoidKernel=Load2Geoid(nmax=self.nmax,deg0scale=0.0)
            self.upliftKernel=Load2Uplift(nmax=self.nmax,deg0scale=0.0)
        else:
            self.geoidKernel=Load2Geoid(knLove=dssnrei.kn,nmax=self.nmax)
            self.upliftKernel=Load2Uplift(hnLove=dssnrei.hn,nmax=self.nmax)

        #setup ocean function
        if p2scache is None:
            p2scache=os.path.join(defaultcache("P2S"),f"p2s_ocean_n{self.nmax}.nc")
        if os.path.exists(p2scache):
            #Read product to sum mat from cache
            shxlogger.info(f"Reading product2sum ocean function from cache: {p2scache}") 
            self.dsp2s_oce=xr.open_dataset(p2scache).cnm.sh.build_nmindex().sh.build_nmindex('_')
        else:
            shxlogger.info(f"Computing ocean function and saving to cache: {p2scache}") 
            if oceansh.sh.nmax != 2*self.nmax:
                oceansh=oceansh.sh.truncate(nmax=self.nmax*2)

            
            self.dsp2s_oce=oceansh.sh.p2s()
            self.dsp2s_oce.sh.drop_nmindex().sh.drop_nmindex('_').to_netcdf(p2scache)
        
        if self.rotfeedback:
            shxlogger.warning("Adding static rotation feedback probably only makes sense for very slowly changing loads (> Chandler wobble freq)")
            self.rotmat=qs_rotfeedback_slow()

    def rotfeed(self,load):
        qsrot=self.rotmat@load
        return qsrot.sh.toggle_nm()

    @staticmethod
    def set_global_mean(load,level):
        load.loc[dict(n=0,m=0)]=level
        return load


    @staticmethod
    def global_mean(load:xr.DataArray):
        """Returns the degree 0, order 0 coefficients of a spherical harmonic dataset"""
        return load.loc[dict(n=0,m=0)]


    def oceanf(self,load=None):
        if load is None:
            #return the ocean function itself
            return self.dsp2s_oce.sel(n_=0,m_=0).drop(['n_','m_','nm_'])
        else:
            #apply the ocean function to a load
            load_oce=self.dsp2s_oce@load
            return load_oce.sh.toggle_nm()

    def load_earth(self,load):
        dsdef=self.geoidKernel(load).to_dataset(name='geoid')
        dsdef['uplift']=self.upliftKernel(load)
        return dsdef
    


