# Test spherical harmonic analysis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025



from fixtures import basin_sim_data 
from shxarray.signal.basinav import Basinav
import numpy as np




def test_sh_basinav(basin_sim_data):
    """ Test (spectral) basin averaging routines"""
    filtername='Gauss500'  
    # filtername='DDK1'  
    #retrieve the basin average on the simulated data
    nmax=int(basin_sim_data.basin_sh.sh.nmax/2)
    datws=basin_sim_data.tws.sh.truncate(nmax)
    engine='shtns' #shtns is faster, but needs to be installed separately
    engine='shlib'
    tws_av=datws.sh.basinav(basin_sim_data.basin_sh,filtername,leakage_corr='scale')
    
    tws_av_no_scale=datws.sh.basinav(basin_sim_data.basin_sh,filtername)
    tws_av_nofilt=datws.sh.basinav(basin_sim_data.basin_sh,leakage_corr='scale')
    tws_av_vishwa2016=datws.sh.basinav(basin_sim_data.basin_sh,filtername,leakage_corr='vishwa2016',engine=engine)
    # import matplotlib.pyplot as plt 
    # ds=basin_sim_data.basin_sh.sum('basins').sh.synthesis(engine=engine)
    # ds.plot()
    # plt.show()
    # breakpoint()
    # fig,axs=plt.subplots(3,1,sharex=True)
    # for ibas in range(3):
        # axs[ibas].plot(basin_sim_data.time,basin_sim_data.basin_avs[:,ibas],label="Truth")
        # axs[ibas].plot(tws_av.time,tws_av[ibas,:],label="scaled")
        # axs[ibas].plot(tws_av_no_scale.time,tws_av_no_scale[ibas,:],label="no scale")
        # axs[ibas].plot(tws_av_vishwa2016.time,tws_av_vishwa2016[ibas,:],label="Vishwakarma2016")
        # axs[ibas].plot(tws_av_nofilt.time,tws_av_nofilt[ibas,:],label="No filter")
    # axs[ibas].legend()
    # plt.show()
    
    #The unfiltered basin average should be close to the truth
    rtol=0.03
    for i in range(basin_sim_data.dims['basins']):
        atol=0.03*np.abs(basin_sim_data.basin_avs[:,i]).max().item()
        closeenough=np.allclose(tws_av_nofilt[i,:],basin_sim_data.basin_avs[:,i],rtol=rtol,atol=atol)
        assert closeenough
    

    #The filtered basin averages can be off by a much larger amount, so we're just going to check if the values are not outrageously off

    rtol=0.5
    for i in range(basin_sim_data.dims['basins']):
        atol=0.5*np.abs(basin_sim_data.basin_avs[:,i]).max().item()
        closeenough=np.allclose(tws_av_no_scale[i,:],basin_sim_data.basin_avs[:,i],rtol=rtol,atol=atol)
        assert closeenough
        closeenough=np.allclose(tws_av[i,:],basin_sim_data.basin_avs[:,i],rtol=rtol,atol=atol)
        assert closeenough

        closeenough=np.allclose(tws_av_vishwa2016[i,:],basin_sim_data.basin_avs[:,i],rtol=rtol,atol=atol)
        assert closeenough




