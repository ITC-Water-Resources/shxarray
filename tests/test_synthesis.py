# Test spherical harmonic synthesis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np
from fixtures import generategrddata,grdvalidation
tol=1e-7

@pytest.fixture
def lon():
    return np.arange(-180,180.01,60.0)


@pytest.fixture
def lat():
    return np.arange(-60,60.01,60.0)

def test_synthesis_full(generategrddata,lon,lat):
    dain,valdata=generategrddata
    dacheck=dain.sh.synthesis(lon,lat)
    #select a subset for checking the input
    dacheck=dacheck.sel(lon=valdata.lon,lat=valdata.lat)
    # # breakpoint()
    assert ((dacheck-valdata)/valdata).max().item() < tol

def test_synthesis_truncated(generategrddata,lon,lat):
    dain,valdata=generategrddata
    dain=dain.sh.truncate(nmin=1)
    dacheck=dain.sh.synthesis(lon,lat)
    #select a subset for checking the input
    dacheck=dacheck.sel(lon=valdata.lon,lat=valdata.lat)
    
    #we can correct the validation data for ignoring n=0, m=0 which has Y00=1 for all positions on the globe
    valdata=valdata-dain.basins*dain.time
    # # breakpoint()
    assert ((dacheck-valdata)/valdata).max().item() < tol


def test_synthesis_1d(generategrddata):
    dain,valdata=generategrddata
    #note: since valdata.lon and valdata.lat share the same dimensions analysis on the pairs is performed rather than on a grid
    dacheck=dain.sh.synthesis(valdata.lon,valdata.lat)
    
    assert ((dacheck-valdata)/valdata).max().item() < tol

def test_synthesis_slice(generategrddata,lon,lat):
    # test synthesis for the case where th input data is a view of the original memory
    dain,valdata=generategrddata
    
    # select a subset from the input and validation
    bassl=slice(1,4)
    timesl=slice(1,3)
    if dain.get_axis_num('nm') == 2:
        dain=dain[bassl,timesl]
    else:
        dain=dain[:,timesl,bassl]

    valdata=valdata[:,timesl,bassl]

    dacheck=dain.sh.synthesis(lon,lat)
    #select a subset for checking the input
    dacheck=dacheck.sel(lon=valdata.lon,lat=valdata.lat)
    assert ((dacheck-valdata)/valdata).max().item() < tol

