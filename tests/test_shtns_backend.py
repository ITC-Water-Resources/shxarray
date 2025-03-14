# Test spherical harmonic and analysis and synthesis of the SHTns backend (https://nschaeff.bitbucket.io/shtns/index.html)
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025


import pytest
import time
import logging
from fixtures import shcapsvalidation,shcapsvalidationaux
import xarray as xr
tol=1e-12

#SHTHNS specific fixture
@pytest.fixture(params=[('C','N'),('C','T'),('F','N'),('F','T')])
def synthinput(shcapsvalidation,request):
    synthdata=shcapsvalidation.sh.synthesis(gtype="gauss",engine="shtns")
   
    #possibly apply request options
    if request.param[1] == 'T':
        synthdata=synthdata.transpose()

    if request.param[0] == 'F':
        synthdata=xr.DataArray(synthdata.data.copy(order='F'),coords=synthdata.coords)
    return synthdata



def test_synthesis_shtns(shcapsvalidation):
    """ Test the SH synthesis-analysis functionality for the SHTns backend (if it is installed),for a Gaussian grid against the SHlib backend"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    engine="shtns"
    #first apply the synthesis on the known spherical caps (using the default grid)
    t0=time.time()
    synthdata=shcapsvalidation.sh.synthesis(engine=engine)
    t1=time.time()
    #also generate a version with shlib for comparison
    synthdata2=shcapsvalidation.sh.synthesis_like(synthdata)

    t2=time.time()
    logging.getLogger().info(f"SHTns synthesis time: {int(1000*(t1-t0))} ms")
    logging.getLogger().info(f"SHlib synthesis time: {int(1000*(t2-t1))} ms")

    dadiff=synthdata-synthdata2
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol)
   
def test_synthesis_shtns_1d(shcapsvalidation):

    """ Test the SH synthesis-analysis functionality for the SHTns backend (for pairs of lonlat) against the SHlib backend"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    engine="shtns"
    #first apply the synthesis on the known spherical caps (using the default grid)
    synthdata=shcapsvalidation.sh.synthesis(lon=shcapsvalidation.plon.data,lat=shcapsvalidation.plat.data,engine=engine)
    #also generate a version with shlib for comparison
    synthdata2=shcapsvalidation.sh.synthesis_like(synthdata)
    
    dadiff=synthdata-synthdata2
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol)
    #also try calling the routine through the synthesis_like functionn
    synthdata3=shcapsvalidation.sh.synthesis_like(synthdata2,engine=engine)
    dadiff=synthdata-synthdata3
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol)



def test_synthesis_shtns_regular(shcapsvalidation):
    """ Test the SH synthesis-analysis functionality for the SHTns backend for a  regular grid against the SHlib backend"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    engine="shtns"
    #first apply the synthesis on the known spherical caps (using the default grid)
    t0=time.time()
    #create a validation on a regular grid starting with lon=0 (compatible with SHTns)
    synthdata2=shcapsvalidation.sh.synthesis(gtype="regular_lon0")

    t1=time.time()
    synthdata=shcapsvalidation.sh.synthesis_like(synthdata2,engine=engine)

    t2=time.time()
    logging.getLogger().info(f"SHLib synthesis time: {int(1000*(t2-t1))} ms")
    logging.getLogger().info(f"SHlib synthesis time: {int(1000*(t1-t0))} ms")
    dadiff=synthdata-synthdata2
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol)

def test_synthesis_shtns_auxdims(shcapsvalidationaux):
    """ Test the SH synthesis-analysis functionality for the SHTns backend for a  regular grid against the SHlib backend. This version has multiple auxiliary dimensions"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    engine="shtns"
    nmax=120 
    #also tests a lower truncation of the SH expansion
    nmin=3
    #first apply the synthesis on the known spherical caps (using the default grid)
    t0=time.time()
    #create a validation on a regular grid starting with lon=0 (compatible with SHTns)
    synthdata2=shcapsvalidationaux.sh.truncate(nmax,nmin).sh.synthesis(engine=engine)

    t1=time.time()
    synthdata=shcapsvalidationaux.sh.truncate(nmax,nmin).sh.synthesis_like(synthdata2)

    t2=time.time()
    logging.getLogger().info(f"SHLib synthesis time: {int(1000*(t2-t1))} ms")
    logging.getLogger().info(f"SHTns synthesis time: {int(1000*(t1-t0))} ms")
    dadiff=synthdata-synthdata2
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol)


def test_analysis_shtns(shcapsvalidation):
    """Test the SHTns analysis functionality"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    
    engine="shtns"
    nmax=100

    synthdata=shcapsvalidation.sh.synthesis(gtype="regular_lon0")
    #now apply the analysis using the SHTnsBackend
    #apply analysis using the SHTns backend
    dsanalysis=synthdata.sh.analysis(nmax=nmax,engine=engine)

    dadiff=shcapsvalidation-dsanalysis

    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<1e-15) 

def test_analysis_shtns_auxdims(shcapsvalidationaux):
    """Test the SHTns analysis functionality"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    
    engine="shtns"
    nmax=100

    synthdata=shcapsvalidationaux.sh.synthesis(gtype="regular_lon0")
    #now apply the analysis using the SHTnsBackend
    #apply analysis using the SHTns backend
    dsanalysis=synthdata.sh.analysis(nmax=nmax,engine=engine)

    dadiff=shcapsvalidationaux-dsanalysis

    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<tol) 

def test_analysis_shtns_order(shcapsvalidation,synthinput):
    """Test the SHTns analysis functionality for various storage orders"""
    try:
        import shtns
        engine="shtns"
    except ImportError:
        pytest.skip("SHTns backend is not installed, skipping test")
    
    engine="shtns"
    nmax=shcapsvalidation.sh.nmax
    #now apply the analysis using the SHTnsBackend
    dsanalysis=synthinput.sh.analysis(nmax=nmax,engine=engine)

    dadiff=shcapsvalidation-dsanalysis

    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff<1e-15) 
