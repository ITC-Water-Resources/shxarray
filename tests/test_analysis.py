# Test spherical harmonic synthesis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import xarray as xr
from shxarray.kernels import ParabolicCap
import os

#high degree testing
nmax=200
tol=1e-7


#low degree testing (quicker)
nmax=20


# validation data set is a bandlimited parabolic cap centered on the positions below
@pytest.fixture
def validation():
    lon=[135 , -71]
    lat=[-31, 65]
    radius=10 #radius of parabolic cap in deg
    pcap=ParabolicCap(nmax,radius)
    pcapsh=pcap.position(lon,lat)
    return pcapsh


# This fixture loads an input grid which is produced by the rltftlbx using
# SH_axisload -ca=10 -l=200 -p=135,-31 | SH_2_nc -I=0.5 -F=shanalysis-test-paracap-n200.nc"
# followed by
# SH_axisload -ca=10 -l=200 -p=-71,65 | SH_2_nc -I=0.5 -F=shanalysis-test-paracap-n200.nc"

@pytest.fixture(params=[('C','N'),('C','T'),('F','N'),('F','T')])
def grdinput(request):
    testdir=os.path.join(os.path.dirname(os.path.realpath(__file__)),'testdata')
    valnc="shanalysis-test-paracap-n200.nc"
    valnc=os.path.join(testdir,valnc)
    dagrd= xr.open_dataset(valnc).z
    
    #possibly apply request options
    if request.param[1] == 'T':
        dagrd=dagrd.transpose()

    if request.param[0] == 'F':
        dagrd=xr.DataArray(dagrd.data.copy(order='F'),coords=dagrd.coords)
    return dagrd

def test_analysis(grdinput,validation):
    """ Test the SH analysis to see if an acceptable sh is computued"""
    dsgrd=grdinput
    checkdata=validation
    dscheck=dsgrd.sh.analysis(nmax)
    # breakpoint()    
    
    #ok we need to rename the input dimension time
    dadiff=checkdata-dscheck.rename(time='nlonlat')
    print("Checking whether retrieved SH solution is within tolerance")
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff < tol)
