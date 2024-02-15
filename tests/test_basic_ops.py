# Test some basic operations of xarray objects filled with spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


import pytest
import shxarray
from shxarray.core.sh_indexing import SHindexBase
import numpy as np

import xarray as xr
from datetime import datetime,timedelta

@pytest.fixture
def sh_sample1():
    startdate=datetime(2023,12,8)
    daterange= [startdate - timedelta(days=dy) for dy in range(60)]
    nmax=30
    auxcoords={"time":daterange}
    da=xr.DataArray.sh.ones(nmax,auxcoords=auxcoords)
    return nmax,0,da


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










