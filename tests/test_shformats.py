# Test some basic operations of xarray objects filled with spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


import pytest
import shxarray
import xarray as xr
import os
import numpy as np
import time
from shxarray.core.logging import shxlogger

def test_icgem_sig():
    """Test ICGEM Backends for xarray"""
    #note gzipped version (should also work)
    icgem1f=os.path.join(os.path.dirname(__file__),'testdata/icgem_test_sig_ITSG.gfc.gz')
    dsicgem=xr.open_dataset(icgem1f,engine='icgem')
    nmax=60
    assert dsicgem.sh.nmax == nmax, "nmax of icgem input is not 60"
    #note some data are missing so we only expect 599 entries
    assert dsicgem.sizes['nm'] == 599, "Not all dat has been read" 
    assert np.datetime64('2004-10-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcnm={(60,46):1.665788238155e-09}
    testsigcnm={(60,-46):2.363195963172e-12}
    for nm,cnm in testcnm.items():
        assert cnm == dsicgem.cnm.loc[:,nm].item(),f"Coefficient {nm} {cnm} not mathching"

    for nm,scnm in testsigcnm.items():
        assert scnm == dsicgem.sigcnm.loc[:,nm].item(),f"Sigma Coefficient {nm} {scnm} not mathching"




def test_icgem_nosig():
    """Test ICGEM Backend for xarray for file without sigmas"""
    icgem2f=os.path.join(os.path.dirname(__file__),'testdata/icgem_test_nosig_ITSG.gfc')
    dsicgem=xr.open_dataset(icgem2f,engine='icgem')
    nmax=2
    assert dsicgem.sh.nmax == nmax, "nmax of icgem input is not 2"
    assert dsicgem.sizes['nm'] == 9, "Not all dat has been read" 
    assert np.datetime64('2013-11-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcnm={(2,0):-4.841696152100e-04,(1,1):0.0}
    for nm,cnm in testcnm.items():
        assert cnm == dsicgem.cnm.loc[:,nm].item(),f"Coefficient {nm} {cnm} not mathching"


def test_icgem_static_tv():
    """Test ICGEM Backend for xarray for file with static and parametric coefficients"""
    icgem3f=os.path.join(os.path.dirname(__file__),'testdata/ITSG-Grace2018s.gfc.gz')
    dsicgem=xr.open_dataset(icgem3f,engine='icgem')
    nmax=200
    assert dsicgem.sh.nmax == nmax, "nmax of icgem input is not 200"
    # assert dsicgem.sizes['nm'] == 9, "Not all dat has been read" 
    # assert np.datetime64('2013-11-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcoef={('cnm_gfct',(10,2)):-9.396504743590e-08,('cnm_gfct',(200,-200)):1.956728673164e-10,
              ('cnm_trnd',(195,159)):7.843896315171e-17,('cnm_trnd',(195,-159)):4.073388066625e-17,
              ('sigcnm_trnd',(195,159)):7.023988334695e-15,('sigcnm_trnd',(195,-159)):7.022127457848e-15,
              ('cnm_acos', (7,3)): 2.269912528926e-11, ('sigcnm_acos',(7,-3)): 4.115793871821e-13,
              ('cnm_acos', (7,-3)): -3.872126870042e-12, ('sigcnm_acos',(7,3)): 4.108227520945e-13}
    # acos    7    3  2.269912528926e-11 -3.872126870042e-12  4.108227520945e-13  4.115793871821e-13 1.0
    for (var,nm),cnm in testcoef.items():
        assert cnm == dsicgem[var].loc[nm,].item(),f"Coefficient {nm} {cnm} not mathching"



