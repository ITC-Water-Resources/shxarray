# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

import xarray as xr
from shxarray.kernels import getSHfilter
from shxarray.signal.leakage_vishwa import leakage_corr_vishwa2016,delta_leakage_corr_vishwa2017

class Basinav:
    def __init__(self, dabasins,filtername=None,leakage_corr=None):

        if leakage_corr is not None and leakage_corr not in ['scale','vishwa2016','vishwa2017']:
            raise RuntimeError(f"Leakage correction method {leakage_corr} not recognized")
        self._leakage_corr=leakage_corr
        self._dabin = dabasins
        self._filtername=filtername

    
    def __call__(self, datws,**kwargs):
        nmax=datws.sh.nmax
        if self._filtername is not None:
            filterOp = getSHfilter(self._filtername,nmax=nmax,transpose=True)
            dabin_f = filterOp(self._dabin.sh.truncate(nmax=nmax))
        else:
            #just truncate to the same nmax as the input
            dabin_f=self._dabin.sh.truncate(nmax)
        
        dabin=self._dabin.sh.truncate(nmax)

        #compute the unscaled basin average
        da_av=(dabin_f@datws)/dabin.sel(n=0,m=0)
        
        if self._leakage_corr in ['scale','vishwa2016']:
            dascales = dabin.dot(dabin,dim='nm')/dabin.dot(dabin_f,dim='nm')
        else:
            dascales=None

        if "engine" in kwargs:
            engine=kwargs["engine"]
        else:
            engine='shlib'
        
        if self._leakage_corr == 'scale':
        
            da_av=da_av*dascales
        elif self._leakage_corr == "vishwa2016":
            leakage=leakage_corr_vishwa2016(datws, dabin, self._filtername,engine=engine) 
            da_av=(da_av-leakage)*dascales
        elif self._leakage_corr == "vishwa2017":
            leakage=leakage_corr_vishwa2016(datws, dabin, self._filtername,engine=engine) 
            deltaleakage=delta_leakage_corr_vishwa2017(datws, dabin, self._filtername,engine=engine)
            #note: no scaling applied
            da_av=(da_av-leakage-deltaleakage)
            # da_av=(da_av-leakage)

        return da_av.drop_vars(['n','m','nm'])

