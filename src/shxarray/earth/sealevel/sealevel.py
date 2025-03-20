# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

import xarray as xr
import numpy as np
from shxarray.core.logging import shxlogger
from shxarray.core.cf import *
from datetime import datetime

class SeaLevelSolver:
    """
    Base Class to aid in solving the self consistent sea level for various discretizations and load types
    """
    def __init__(self,rotfeedback=False):
        #whether to apply a simple (static) rotational feedback
        self.rotfeedback=rotfeedback
    
    @staticmethod
    def set_global_mean(load,level):
        """Sets the global mean value of of a load. Depending on the discretization, this needs to be differently implemented"""
        
        raise NotImplementedError("set_global_mean(), needs to be implemented in derived class")
    
    @staticmethod
    def global_mean(load):
        raise NotImplementedError("global_mean(), needs to be implemented in derived class")
    def rotfeed(self,load):
        raise NotImplementedError("rotfeed(), needs to be implemented in derived class")

    def oceanf(self,load=None):
        raise NotImplementedError("oceanf(), needs to be implemented in derived class")

    def load_earth(self,load):
        raise NotImplementedError("load_earth(), needs to be implemented in derived class")
    
    

    def __call__(self, load_force):
        """Iteratively solve the sea level equation for a given load"""  
        
        #compute the mean contributions of the fixed (forced) load
        force00=self.global_mean(load_force)
        #Set up initial ocean load ocean function
        unioce=self.oceanf()
        oce_surf_ratio=self.global_mean(unioce)
        dphi_g=-force00/oce_surf_ratio
        #initial value for the ocean load
        load_sea=dphi_g*unioce
        
        relratio=1
        rel_thres=1e-5
        maxit=7
        it=0
        damp=0.95 
        while it < maxit and rel_thres < abs(relratio):
            #compute the loading response
            load_tot=load_sea+load_force
            ds_loadresp=self.load_earth(load_tot)
            quasi_sea=ds_loadresp.geoid-ds_loadresp.uplift
            if self.rotfeedback:
                qsrot=self.rotfeed(load_tot)
                #add static rotational feedback component to quasi_sea
                quasi_sea.loc[qsrot.nm]+=qsrot
            self.set_global_mean(quasi_sea,dphi_g) 
            #compute new ocean load
            load_sea=self.oceanf(quasi_sea)
            delta=force00+self.global_mean(load_sea)
            relratio=np.max(np.abs(delta/force00)).item()
            #update dphi_g in the right direction
            dphi_g-=damp*delta/oce_surf_ratio
            shxlogger.info(f"current mass inconsistency: iteration {it}, relratio:{relratio}")
            
            it+=1

        
        quasi_sea.name="quasi_sea"
        load_sea.name="load_sea"
        load_force.name="load_force"
        ds_loadresp.geoid.name="geoid"
        ds_loadresp.uplift.name="uplift"

        dsout=xr.merge([quasi_sea,load_sea,load_force,ds_loadresp.geoid,ds_loadresp.uplift])
        
        #add some useful attributes to the output
        dsout.attrs['history'] = str(datetime.utcnow()) + f": shxarray ({self.__class__.__name__}"
        

        dsout.quasi_sea.attrs["long_name"]="Quasi spectral sea level (non zero over land)"
        dsout.quasi_sea.attrs["units"]="m"
        
        dsout.load_sea.attrs["long_name"]="Relative Sea level (load) as equivalent water height"
        dsout.load_sea.attrs["units"]="m"
        
        dsout.load_force.attrs["long_name"]="Applied load as equivalent water height"
        dsout.load_force.attrs["units"]="m"

        dsout.geoid.attrs["long_name"]="Geoid height change induced by combined load (force+sea)"
        dsout.geoid.attrs["units"]="m"

        dsout.uplift.attrs["long_name"]="Uplift induced by combined load (force+sea)"
        dsout.uplift.attrs["units"]="m"


        return dsout
        


        
