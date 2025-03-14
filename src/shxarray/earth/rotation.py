# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
# Source: See chapter 2.4.3 of Rietbroek et al 2014 http://nbn-resolving.de/urn:nbn:de:hbz:5n-35460i and corrections in https://github.com/strawpants/rlftlbx/blob/master/SHtlbx/rotfeedback.f90

import xarray as xr
import numpy as np
from shxarray.earth.constants import a_earth, rho_water, ohm_chandler,ohm_earth,ixx_A,izz_C,k2loadprem, k2bodyprem,h2bodyprem, g

from shxarray.core.sh_indexing import SHindexBase

def chi_T2J():
    fac=np.pi*(a_earth**4)*rho_water
    chimat=np.zeros([3,4])
    chimat[0,2]=-4.0/5.0*np.sqrt(10.0/6.0)
    chimat[1,3]=chimat[0,2]
    chimat[2,0]=8.0/3.0
    chimat[2,1]=-8/(3*np.sqrt(15))
    return fac*chimat
    
def gamma_j2m(k2load=None):
    if k2load is None:
        k2load=k2loadprem

    gmat=np.zeros([3,3])
    gmat[0,0]=(ohm_earth/ohm_chandler)*(1+k2load)/ixx_A
    gmat[1,1]=gmat[0,0]
    gmat[2,2]=-(1+k2load)/izz_C

    return gmat


def psi_m2ilam():
    pmat=np.zeros([4,3])
    pmat[0,2]=2/3
    pmat[1,2]=-2/(3*np.sqrt(5))
    pmat[2,0]=-1/(np.sqrt(15))
    pmat[3,1]=pmat[2,0]
    fac=(a_earth*ohm_earth)**2
    return fac*pmat

def t_lam2s():
    tmat=np.zeros([3,4])
    tmat[0,1]=1
    tmat[1,2]=1
    tmat[2,3]=1

    fac=(1+k2bodyprem-h2bodyprem)/g
    return fac*tmat

def qs_rotfeedback_slow():
    """
    Returns a matrix relating sh surface load coefficients to degree 2 changes in Quasi spectral sea level
    """

    chi=chi_T2J()
    gam=gamma_j2m()
    psi=psi_m2ilam()
    tm=t_lam2s()
    rotmat=tm@psi@gam@chi
    nmi_=SHindexBase.mi_fromtuples([(2,0),(2,1),(2,-1)],'_')
    nmi=SHindexBase.mi_fromtuples([(0,0),(2,0),(2,1),(2,-1)])
    return xr.DataArray(rotmat,dims=["nm_","nm"],coords=[nmi_,nmi]) 


    





