# Test Wigner 3j symbols and Gaunt  
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np
from sympy.physics.wigner import wigner_3j,gaunt,real_gaunt
import os

tol=1e-16


# def get_sympyvaldata():
    # """Create a on-disk dataset with validation coefficients for Wigner3j, Gaunt and realGaunt coefficients"""
    # sympyfile=os.path.join(os.path.dirname(__file__),'testdata/sympy_valdata.pkl')
    # sympyval={}
    # gauntreal={}
    # precrg=3
    # if not os.path.exists(sympyfile):
        # #generate new data
        # nmax=4
        # for n2 in range(0,nmax+1):
            # for n3 in range(0,nmax+1):
                # for m2 in range(0,n2+1):
                    # for m3 in range(0,n3+1):
                        # n1minrg=max(abs(n2-n3),abs(m2+m3))
                        # n1minrg+=(n2+n3+n1minrg)%2
                        # n1maxrg=n2+n3

                        # # gauntreal[(n2,n3,m2,m3)]=[float(val) for val in [real_gaunt(n1,n2,n3,m1,m2,m3,prec=precrg) for n1 in n]]
                          




def generate_w3j_validation(jmin,jmax,j2,j3,m1,m2,m3):
    #note this requires the python module sympy
    w3jval=[float(val) for val in [wigner_3j(j1,j2,j3,m1,m2,m3) for j1 in range(jmin,jmax+1)]]
    return w3jval

def test_Wigner3j():
    for j2 in [16,41,300]:
        for j3 in [13,121]:
            for m2 in [-2,100]:
                for m3 in [13,5]:
                    if j2 < m2 or j3 < m3:
                        continue

                    daw3j=xr.DataArray.sh.wigner3j(j2,j3,m2,m3)
                    #note only one valid m1 should be present 

                    m1=daw3j.m.data[0]
                    w3jval=generate_w3j_validation(daw3j.j.min().item(),daw3j.j.max().item(),j2,j3,m1,m2,m3)
                    assert(np.allclose(daw3j.data,w3jval,rtol=tol))


def generate_gaunt_validation(n,n2,n3,m1,m2,m3):
    #note this requires the python module sympy
    gauntval=[float(val) for val in [gaunt(n1,n2,n3,m1,m2,m3) for n1 in n]]
    return gauntval

def test_Gauntnormal():
    # for n2 in [10,11,300]:
        # for n3 in [90,121]:
            # for m2 in [-11,100]:
                # for m3 in [31,5]:
    nmax=8
    for n2 in range(0,nmax+1):
        for n3 in range(0,nmax+1):
            for m2 in range(0,n2+1):
                for m3 in range(0,n3+1):
                    if n2 < m2 or n3 < m3:
                        continue
                    dagaunt=xr.DataArray.sh.gaunt(n2,n3,m2,m3)
                    m1=dagaunt.m.data[0]

                    gauntval=generate_gaunt_validation(dagaunt.n.data,n2,n3,m1,m2,m3)
                    assert(np.allclose(dagaunt.data,gauntval,rtol=tol))

def generate_gauntreal_validation(n,n2,n3,m1,m2,m3):
    #note this requires the python module sympy
    prec=3
    gauntval=[float(val) for val in [real_gaunt(n1,n2,n3,m1,m2,m3,prec=prec) for n1 in n]]
    return gauntval

def test_Gauntreal():
    nmax=5
    for n2 in range(0,nmax+1):
        for n3 in range(0,nmax+1):
            for m2 in range(0,n2+1):
                for m3 in range(0,n3+1):
                    if n2 < m2 or n3 < m3:
                        assert(False)
                        continue
                
                    dagaunt=xr.DataArray.sh.gauntReal(n2,n3,m2,m3)
                    ms=np.unique(dagaunt.m)
                    gauntval=[]
                    # breakpoint()
                    for m1 in ms:
                        gauntval.extend(generate_gauntreal_validation(dagaunt.sel(m=m1).n.data,n2,n3,m1,m2,m3))
                    
                    # breakpoint()
                    closeEnough=np.allclose(dagaunt.data,gauntval,rtol=1e-3)
                    if not closeEnough:
                        breakpoint()
                    assert(closeEnough)
                    
