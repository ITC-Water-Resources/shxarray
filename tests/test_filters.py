# Test DDK filtering of spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024



import pytest
import xarray as xr
import os
import numpy as np


#relative error
tol=1e-11


# #low degree testing (quicker)
# nmax=20


# Validation data set ( see https://github.com/strawpants/GRACE-filter/tree/master/tests) 
@pytest.fixture
def shinput():
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    insh=os.path.join(testdir,"GSM-2_2008122-2008153_0030_EIGEN_G---_0004in.gz")
    ds=xr.open_dataset(insh,engine='shascii')
    return ds

@pytest.fixture
def shoutput():
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    outsh=os.path.join(testdir,"GSM-2_2008122-2008153_0030_EIGEN_G---_0004out.gz")
    ds=xr.open_dataset(outsh,engine='shascii')
    return ds

@pytest.fixture
def shoutputn60():
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    outsh=os.path.join(testdir,"GSM-2_2008122-2008153_0030_EIGEN_G---_0004lmax60out.gz")
    ds=xr.open_dataset(outsh,engine='shascii')
    return ds




def test_ddkbasic(shinput,shoutput):
    """ Test the filtering of SH coefficients using an anisotropic ddk filter"""
    dacheck=shinput.cnm.sh.filter('DDK2')
    dadiff=(dacheck-shoutput.cnm)/shoutput.cnm
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff < tol)
    

def test_ddknmax60(shinput,shoutputn60):
    """ Test the filtering of  a truncated set of SH coefficients using an anisotropic ddk filter"""
    shin=shinput.cnm.sh.truncate(60)
    dacheck=shin.sh.filter('DDK2')
    dadiff=(dacheck-shoutputn60.cnm)/shoutputn60.cnm
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff < tol)

def test_ddkmult(shinput,shoutput):
    """Test the filtering on a set (multiple) of SH coefficients using an anisotropic ddk filter"""
    
    nbasins=5    

    basinscales=np.arange(nbasins)

    scales1=xr.DataArray(basinscales,coords={"basins":basinscales})
    shinm=shinput.cnm*scales1

    shoutm=shoutput.cnm*scales1
    
    dacheck=shinm.sh.filter('DDK2')
    dadiff=(dacheck-shoutm)/shoutm
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff < tol)

@pytest.fixture
def shgausstest():
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    outsh=os.path.join(testdir,"gauss300testoutsub.sh.gz")
    ds=xr.open_dataset(outsh,engine='shascii')
    return ds

def test_Gauss(shinput,shgausstest):
    shfilt=shinput.cnm.sh.filter("Gauss300")
    
    dadiff=(shgausstest.cnm-shfilt.loc[shgausstest.nm])/shgausstest.cnm
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    reltol=1e-7
    assert(maxdiff < reltol)


