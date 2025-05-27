# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

import xarray as xr
import numpy as np
from scipy.signal import hilbert
from shxarray.kernels import getSHfilter
from shxarray.core.logging import shxlogger

def leakage_corr_vishwa2016(datws, dabasins, filtername,engine='shlib'):
    """
    Basin averages leakage correction as described in Vishwakarma et al. (2016) for the given Total water storage dataset and basins.
    source: https://onlinelibrary.wiley.com/doi/abs/10.1002/2016WR018960

    #note this only provides the Ic_hat correction (in eq. 36) not the estimate of the basin average itself

    Parameters
    ----------
    datws : xarray.Dataset
        Dataset containing the terrestrial water storage data.
    dabasins : xarray.DataArray
        Basin representation of the data
    filtername : str
        Name of the smoothing filter to apply
    
    Returns
    -------
    xarray.Dataset
        Dataset containing basin averaged leakage correction data.
    """
    timedim='time'
    #check if a time dimension is present int he input data
    if timedim not in datws.dims:
        raise RuntimeError("No time dimension found in the input data")

    nbasindims=len(dabasins.dims)
    if nbasindims == 1:
        dabasins=dabasins.expand_dims('basins')
    elif nbasindims != 2:
        raise RuntimeError("Basins dimensions can only be 1 or 2")

    #compute the complementary basins
    basins_comp=-dabasins
    basins_comp.loc[dict(n=0,m=0)]+=1
    
    filterOp_t = getSHfilter(filtername,nmax=dabasins.sh.nmax,transpose=True)
    #compute Knm coefficients (note that we use the transpose of the filter, when relevant)
    daknm=basins_comp.sh.multiply(filterOp_t(dabasins),engine=engine,truncate=True)
    
    #equidistant sampling interval and make sure it's even
    delta_t=datws[timedim].diff(dim=timedim).median().item()
    ntime=int((datws[timedim].max().item()-datws[timedim].min().item())/delta_t)
    #mkae sure it's even (for fft purposes)
    ntime=ntime+ntime%2
    time_eq=np.linspace(datws[timedim].min().item(),datws[timedim].max().item(),ntime,dtype=datws[timedim].dtype)

    #compute single filtered leakage signal (and sample at equidistant time)
    filterOp = getSHfilter(filtername,nmax=datws.sh.nmax)
    datws_f=filterOp(datws)
    
     
    Ic_f=((daknm@datws_f)/dabasins.sel(n=0,m=0)).interp({timedim:time_eq})
    #compute double filtered leakage signal
    datws_ff=filterOp(datws_f)
    Ic_ff=((daknm@datws_ff)/dabasins.sel(n=0,m=0)).interp({timedim:time_eq})
    
    #compute the phase shift between Ic_f and Ic_ff using Phillips et al 2012 (e.q. 4, 6)
    #i.e. solve with least squares a,b,c:  Ic_f = a + b * Ic_ff + c * imag(Hilbert(Ic_ff))
    
    
    auxdim=[d for d in Ic_ff.dims if d != 'time'][0]
    #allocate design matrix
    A=np.zeros([len(time_eq),3])

    phase_exp=[]
    for i in range(dabasins.sizes[auxdim]):
        A[:,0]=1
        A[:,1]=Ic_ff.isel({auxdim:i}).data
        A[:,2]=np.imag(hilbert(Ic_ff.isel({auxdim:i}).data))
        try:
            reg_fit,_,_,_=np.linalg.lstsq(A,Ic_f.isel({auxdim:i}).data)
        except np.linalg.LinAlgError:
            shxlogger.warning(f"Least squares fit failed for basin {i}, phase for leakage could not be computed, assuming np phase shift")
            reg_fit=[0,1,0]
        # compute phase=atan(c/b) to get the phase shift and take the complex exponential
        phase_exp.append(np.exp(-1j*(np.arctan2(reg_fit[2],reg_fit[1]))))
    
    
    ax=Ic_ff.dims.index('time')
    #shift Ic_ff towards Ic_f using fft -> ifft
    fftdat=(phase_exp*np.fft.rfft(Ic_ff.data,axis=ax).T).T
    Ic_ff_shft=np.fft.irfft(fftdat,axis=ax)
    
    
    #apply the correction to the single filtered series
    fftdat=(phase_exp*np.fft.rfft(Ic_f.data,axis=ax).T).T
    Ic_f_shft=xr.zeros_like(Ic_f)
    Ic_f_shft[:,:]=np.fft.irfft(fftdat,axis=ax)
    #compute the median ratio between the shifted double filtered series and the single filtered series
    Ic_frac=(Ic_f/Ic_ff_shft).mean(timedim)
    
    
    #resample to original time
    Ic_hat=Ic_f_shft.interp({timedim:datws[timedim]})
    #rescale
    Ic_hat=Ic_hat*Ic_frac
    return Ic_hat

# def delta_kernels_vishwa2017(dabasins):
    # """
    # Compute the kernels which can be used for computing the delta compo 
    # """

def delta_leakage_corr_vishwa2017(datws, dabasins, filtername,engine='shlib'):
    """
    Provides the delta leakage correction as described in Vishwakarma et al. (2017) for the given Total water storage dataset and basins.
    source: https://onlinelibrary.wiley.com/doi/abs/10.1002/2017WR021150

    #note this only provides an estimate of the deltaFc correction part (in eq. 13) not the estimate of the basin average itself
    
    deltaF is defined as deltaF = F - fc with fc the true catchment average and F the masked true field within the basin.
    The masking operation can be written as 
    F = B f and its complement F* = (I-B) f, where B is the basin operator and f the true field.

    When averaging deltaFwith the basin averaging operator b' (basin vectors divided by their area) we can rewrite the above as a dedicated basin defined delta operator
    
    deltaF_ave = b' B f - b' f =  deltaOp f 
    with deltaOp = b' (B-I)
    In the spectral deomain deltaOp is equivalent to a transposed spherical harmonic vector, independent on the field f and the filter operation 
    
    When applying once (_f) and twice (_ff) filtered versions we can thus compute
    deltaF_ave_f= deltaOp f_f
    and 
    deltaF_ave_ff= deltaOp f_ff

    An approximation of the true deltaF_ave can be computed as

    deltaF_ave = ratio* deltaF_ave_f 

    where the ratio comes from the median values of the once/twice filtered deltaF_ave's
    ratio =median(deltaF_ave_f(ti)/deltaF_ave_ff(ti)) with ti spanning the times


    Parameters
    ----------
    datws : xarray.Dataset
        Dataset containing the terrestrial water storage data.
    dabasins : xarray.DataArray
        Basin representation of the data
    filtername : str
        Name of the smoothing filter to apply
    
    Returns
    -------
    xarray.Dataset
        Dataset containing basin averaged leakage correction data.
    """
    
    #filter data once/twice

    filterOp = getSHfilter(filtername,nmax=datws.sh.nmax)
    datws_f=filterOp(datws)
    datws_ff=filterOp(datws_f)

    nmax=datws.sh.nmax
    #Complementary basins (representing the outside of the catchments)
    #compute the complementary basins
    basins_comp=-dabasins.sh.truncate(nmax)
    basins_comp.loc[dict(n=0,m=0)]+=1
    
    #mean
    Ac=dabasins.sel(n=0,m=0)

    #delta operator deltaOp= b' (B-I) = ((B-I)' b)'  
    deltaOp=basins_comp.sh.multiply(dabasins.sh.truncate(nmax),engine=engine,truncate=True)/Ac
    breakpoint()
    #apply the delta operator to the once/twice filtered data
    delta_ave_f= deltaOp@datws_f
    delta_ave_ff= deltaOp@datws_ff

    #retrieve the median ratio over the time dimension (assumes 'time' is a valid dimension in datws)
    deltaratio=(delta_ave_f/delta_ave_ff).median('time')
    
    #return appximate deltaF_ave
    return delta_ave_f*deltaratio





