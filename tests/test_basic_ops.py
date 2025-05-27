# Test some basic operations of xarray objects filled with spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


import pytest
import shxarray
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.logging import shxlogger
import numpy as np

import xarray as xr
from datetime import datetime,timedelta
import time

@pytest.fixture
def sh_sample1():
    startdate=datetime(2023,12,8)
    daterange= [startdate - timedelta(days=dy) for dy in range(60)]
    nmax=30
    auxcoords={"time":daterange}
    da=xr.DataArray.sh.ones(nmax,auxcoords=auxcoords)
    #make coefficients decline with degree
    power=2
    degscale=xr.DataArray([1/((n+1)**power) for n in da.n],dims="nm")
    return nmax,0,da*degscale


@pytest.fixture
def sh_truncated(sh_sample1):
    _,_,dasample=sh_sample1
    nmax2=28
    nmin2=3
    datrunc=dasample.sh.truncate(nmax2,nmin2)
    return nmax2,nmin2,datrunc

def test_truncate(sh_truncated):
    nmax2,nmin2,datrunc=sh_truncated
    assert datrunc.sh.nmax == nmax2, "nmax of truncated DataArray does not agree with expectation"
    assert datrunc.sh.nmin == nmin2, "nmin of truncated DataArray does not agree with expectation"
    assert len(datrunc.nm) == SHindexBase.nsh(nmax2,nmin2,squeeze=True), "Size of truncated DataArray does not agree with expectation"


def test_add_sub(sh_sample1,sh_truncated):
    nmax,nmin,dasample=sh_sample1
    nmax2,nmin2,datrunc=sh_truncated
    da_add=datrunc + 2*dasample
    assert da_add.sh.nmin == nmin2, "Nmin of adding results is not consistent"
    assert da_add.sh.nmax == nmax2, "Nmax of adding results is not consistent"
    assert da_add.sum().item() == 3*np.prod(da_add.shape), "Adding operation not consistent"

    da_sub=dasample - 2* datrunc
    assert da_sub.sh.nmin == nmin2, "Nmin of subtracting results is not consistent"
    assert da_sub.sh.nmax == nmax2, "Nmax of subtracting results is not consistent"
    assert da_sub.sum().item() == -np.prod(da_sub.shape), "Adding operation not consistent"


@pytest.mark.skip(reason="Not yet properly implemented")
def test_multiply(sh_sample1,sh_truncated):
    
    nmax1,nmin,dasample=sh_sample1
    nmax2,nmin2,datrunc=sh_truncated
    # #truncate
    # dasample=dasample.sh.truncate(nmax1)
    # datrunc=datrunc.sh.truncate(nmax2)
    
    #multiply by 2
    datrunc*=2
    # t0=time.time()
    # daout=dasample.sh.multiply(datrunc,engine="shtns")
    t1=time.time()

    # shxlogger.info(f"Time to multiply with shtns {t1-t0:.2f} seconds")
    # breakpoint()
    daout2=dasample.sh.multiply(datrunc,engine="shlib",method="spatial")
    t2=time.time()
    shxlogger.info(f"Time to multiply with shlib spatial {t2-t1:.2f} seconds")

    # daout=dasample.sh.multiply(datrunc,engine="exp")
    # breakpoint()

    daout3=dasample.sh.multiply(datrunc,engine="shlib",method="spectral")
    t3=time.time()
    shxlogger.info(f"Time to multiply with shlib (spectral) {t3-t2:.2f} seconds")
    
    # dsout=daout.to_dataset(name='shtns')
    # dsout['shlibspat']=daout2
    # dsout['shlibspec']=daout3
    # dsout.sh.synthesis().to_netcdf('tests/testdata/mult_test.nc',mode='w')
    

    assert np.allclose(daout,daout2.sh.truncate(nmax2),atol=1e-13), "Results of multiply with different engines do not agree" 







