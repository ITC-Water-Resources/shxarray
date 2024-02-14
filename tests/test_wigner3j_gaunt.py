# Test Wigner 3j symbols and Gaunt  
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np
import gzip
import os
import pickle


@pytest.fixture
def wigner3jVal():
    """Create a on-disk dataset with validation coefficients for Wigner3j symbols"""
    sympyfile=os.path.join(os.path.dirname(__file__),'testdata/sympy_wigner3jvalidation.pkl.gz')

    wigner3j=[]
    if not os.path.exists(sympyfile):
        from sympy.physics.wigner import wigner_3j
        #generate new data
        nmax=9
        for n2 in range(0,nmax+1):
            for n3 in range(0,nmax+1):
                for m2 in range(-n2,n2+1):
                    for m3 in range(-n3,n3+1):
                        #retrieve non-zero orders and degrees
                        daw3j=xr.DataArray.sh.wigner3j(n2,n3,m2,m3)
                        nm=daw3j.jm.data
                        m1=np.unique(daw3j.m)[0]
                        data=[float(val) for val in [wigner_3j(n1,n2,n3,m1,m2,m3) for n1 in daw3j.j]]
                        wigner3j.append({"nm":nm.copy(),"n2":n2,"n3":n3,"m2":m2,"m3":m3,"data":data.copy()})
                        # breakpoint()
        with gzip.open(sympyfile,'wb') as fid:
            pickle.dump(wigner3j,fid,protocol=1)
    else:
        with gzip.open(sympyfile,'rb') as fid:
            wigner3j=pickle.load(fid)
    
    return wigner3j



def test_W3j(wigner3jVal):
    rtol=1e-8
    for valdata in wigner3jVal:
        daw3j=xr.DataArray.sh.wigner3j(valdata["n2"],valdata["n3"],valdata["m2"],valdata["m3"])
        assert(len(valdata["data"]) == len(daw3j))
        closeEnough=np.allclose(daw3j.data,valdata["data"],rtol=rtol)
        assert(closeEnough)



@pytest.fixture
def gauntVal():
    """Create a on-disk dataset with validation coefficients for Gaunt coefficients"""
    sympyfile=os.path.join(os.path.dirname(__file__),'testdata/sympy_gauntvalidation.pkl.gz')

    gauntv=[]
    if not os.path.exists(sympyfile):
        from sympy.physics.wigner import gaunt
        #generate new data
        nmax=9
        for n2 in range(0,nmax+1):
            for n3 in range(0,nmax+1):
                for m2 in range(-n2,n2+1):
                    for m3 in range(-n3,n3+1):
                        #retrieve non-zero orders and degrees
                        dagaunt=xr.DataArray.sh.gaunt(n2,n3,m2,m3)
                        nm=dagaunt.nm.data
                        m1=np.unique(dagaunt.m)[0]
                        data=[float(val) for val in [gaunt(n1,n2,n3,m1,m2,m3) for n1 in dagaunt.n.data]]
                        gauntv.append({"nm":nm.copy(),"n2":n2,"n3":n3,"m2":m2,"m3":m3,"data":data.copy()})
        with gzip.open(sympyfile,'wb') as fid:
            pickle.dump(gauntv,fid,protocol=1)
    else:
        with gzip.open(sympyfile,'rb') as fid:
            gauntv=pickle.load(fid)
    
    return gauntv


def test_Gauntnormal(gauntVal):
    rtol=1e-8
    for valdata in gauntVal:
        dagaunt=xr.DataArray.sh.gaunt(valdata["n2"],valdata["n3"],valdata["m2"],valdata["m3"])
        assert(len(valdata["data"]) == len(dagaunt))
        closeEnough=np.allclose(dagaunt.data,valdata["data"],rtol=rtol)
        assert(closeEnough)


@pytest.fixture
def gauntrealVal():
    """Create a on-disk dataset with validation coefficients for real Gaunt coefficients"""
    sympyfile=os.path.join(os.path.dirname(__file__),'testdata/sympy_realgauntvalidation.pkl.gz')

    gauntv=[]
    if not os.path.exists(sympyfile):
        from sympy.physics.wigner import real_gaunt
        #generate new data
        nmax=8
        for n2 in range(0,nmax+1):
            for n3 in range(0,nmax+1):
                for m2 in range(-n2,n2+1):
                    for m3 in range(-n3,n3+1):
                        #retrieve non-zero orders and degrees
                        dagaunt=xr.DataArray.sh.gauntReal(n2,n3,m2,m3)
                        if len(dagaunt) == 0:
                            # no-nonzero values
                            continue
                        nm=dagaunt.nm.data
                    
                        m1=np.unique(dagaunt.m)[0]
                        data=[float(val) for val in [real_gaunt(n1,n2,n3,m1,m2,m3) for n1 in dagaunt.n.data]]
                        gauntv.append({"nm":nm.copy(),"n2":n2,"n3":n3,"m2":m2,"m3":m3,"data":data.copy()})
        with gzip.open(sympyfile,'wb') as fid:
            pickle.dump(gauntv,fid,protocol=1)
    else:
        with gzip.open(sympyfile,'rb') as fid:
            gauntv=pickle.load(fid)
    
    return gauntv


def test_Gauntreal(gauntrealVal):
    rtol=1e-8
    for valdata in gauntrealVal:
        dagaunt=xr.DataArray.sh.gauntReal(valdata["n2"],valdata["n3"],valdata["m2"],valdata["m3"])
        assert(len(valdata["data"]) == len(dagaunt))
        closeEnough=np.allclose(dagaunt.data,valdata["data"],rtol=rtol)
        # m1=dagaunt.m.data[0]
        # m2=valdata["m2"]
        # m3=valdata["m3"]
        # if not closeEnough:
            # print(f'not ok {np.sign(m1*m2*m3)}')            
        # else:
            # print(f'OK {np.sign(m1*m2*m3)}')            
        assert(closeEnough)
