# Test comversion of gravity functionals in the spherical harmonic domain
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025



import pytest
import xarray as xr
import os
import numpy as np


#relative error
tol=1e-11

@pytest.fixture
def shinput():
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    insh=os.path.join(testdir,"GSM-2_2008122-2008153_0030_EIGEN_G---_0004in.gz")
    ds=xr.open_dataset(insh,engine='shascii')
    return ds


def test_stokes_geoid_tws_stokes(shinput):
    """Self consistency test for gravity functionals stokes -> geoid -> tws -> stokes"""
    
    try:
        dageoid=shinput.cnm.sh.geoid()
        assert False, "geoid conversion should not be available when input type is not defined"
    except TypeError:
        #ok as the conversion should fail when the input gravtype is not set
        pass

    #either use 
    dageoid=shinput.cnm.sh.geoid(ingravtype="stokes")
    #or set the gravtype on the input
    shinput.cnm.sh.gravtype='stokes'
    dageoid=shinput.cnm.sh.geoid()
    
    datws=dageoid.sh.tws()
    
    #convert back to stokes
    dastokes=datws.sh.stokes()

    #we should actually end up with what we had before
    assert np.allclose(shinput.cnm,dastokes)


def test_stokes_tws_uplift_stokes(shinput):
    """Self consistency test for gravity functionals stokes -> tws -> uplift -> stokes"""

    

    datws=shinput.cnm.sh.tws(ingravtype="stokes")
    dauplift=datws.sh.uplift()  
    
    #convert back to stokes
    dastokes=dauplift.sh.stokes()

    #we should actually end up with what we had before
    assert np.allclose(shinput.cnm,dastokes)

    
def test_stokes_horzdef_geoid_stokes(shinput):
    """Self consistency test for gravity functionals stokes -> horzdef -> geoid -> stokes"""
    
    shinput.cnm.sh.gravtype='stokes'
    dahorzdef=shinput.cnm.sh.horzdef()
    dageoid=dahorzdef.sh.geoid()  
    
    #convert back to stokes
    dastokes=dageoid.sh.stokes()

    #we should actually end up with what we had before
    assert np.allclose(shinput.cnm,dastokes)
