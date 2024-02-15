# Test spherical harmonic synthesis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np

nmax=20
tol=1e-7

#high degree testing
nmax=200


# The validation datatset is produced from the RLFTLBX using the following commands:
# SH_dump -l200 -u | SH_2_nc  -R=-180/180/-60/60 -I=60 -g
# Which creates a unit spherical harmonic dataset mapped to a global grid
@pytest.fixture
def validation():
    # holds for nmax=200
    if nmax == 200:
    #validation data constructed with rlfxtblx tools: SH_dump -l200 -u | SH_2_nc  -R=-180/180/-60/60 -I=60 -g
        valdata=[(-180,60,-2.48350572586),(-120,60,-2.68359804153),(-60,60,22.2082977295),(0,60,3538.09448242),(60,60,95.0010375977),(120,60,-6.011282444),(180,60,-2.48350572586),(-180,0,2.33900880814),(-120,0,-3.7789080143),(-60,0,-7.44546556473),(0,0,559.570007324),(60,0,3.45617771149),(120,0,-0.587444782257),(180,0,2.33900880814),(-180,-60,12.755859375),(-120,-60,6.24348688126),(-60,-60,-0.0657368898392),(0,-60,0.941503584385),(60,-60,1.90674889088),(120,-60,-21.6708869934),(180,-60,12.755859375)]
    elif nmax == 20:
    #validation data constructed with rlfxtblx tools: SH_dump -l20 -u | SH_2_nc  -R=-180/180/-60/60 -I=60 -g
        valdata=[(-180,60,-2.52639126778),(-120,60,-3.54197382927),(-60,60,-26.5734081268),
                 (0,60,121.698486328),(60,60,26.0581665039),(120,60,-5.30710315704),
                 (180,60,-2.52639126778),(-180,0,1.5869371891),(-120,0,-2.11619329453),
                 (-60,0,-4.59407281876),(0,0,33.3332214355),(60,0,1.56787621975),
                 (120,0,0.476241528988),(180,0,1.5869371891),(-180,-60,3.8463101387),
                 (-120,-60,-4.64725017548),(-60,-60,-0.319970816374),(0,-60,0.918827295303),
                 (60,-60,1.72289800644),(120,-60,-6.23128652573),(180,-60,3.8463101387)]
    lon=[lon for lon,_,_ in valdata]
    lat=[lat for _,lat,_ in valdata]
    values=[val for _,_,val in valdata]

    valdata=xr.DataArray(values,coords={"lon":("nlonlat",lon),"lat":("nlonlat",lat)},dims=["nlonlat"])
    

    return valdata


# this fixture produces a spherical unit dataset of with multiple dimensions in both c and f storage order variants
@pytest.fixture(params=[("C","N"),("C","T"),("F","N"),("F","T")])
def generatedata(validation,request):
    nbasins=5
    ntime=4
    basinscales=np.arange(nbasins)
    timescales=0.3*np.arange(ntime)

    scales1=xr.DataArray(basinscales,coords={"basins":basinscales})
    scales2=xr.DataArray(timescales,coords={"time":timescales})
    coords={"basins":basinscales,"time":timescales}
    if request.param[0] == 'C': 
        dain=xr.DataArray.sh.ones(nmax ,auxcoords=coords,order='C')*scales1*scales2
    else:
        dain=xr.DataArray.sh.ones(nmax ,auxcoords=coords,order='F')*scales1*scales2
    
    if request.param[1] == 'T':
        dain=dain.transpose()
    

    #also scale the validation data
    valdata=validation*scales1*scales2

    return dain,valdata

@pytest.fixture
def lon():
    return np.arange(-180,180.01,60.0)


@pytest.fixture
def lat():
    return np.arange(-60,60.01,60.0)

def test_synthesis_full(generatedata,lon,lat):
    dain,valdata=generatedata
    dacheck=dain.sh.synthesis(lon,lat)
    #select a subset for checking the input
    dacheck=dacheck.sel(lon=valdata.lon,lat=valdata.lat)
    # # breakpoint()
    assert ((dacheck-valdata)/valdata).max().item() < tol

def test_synthesis_truncated(generatedata,lon,lat):
    dain,valdata=generatedata
    dain=dain.sh.truncate(nmin=1)
    dacheck=dain.sh.synthesis(lon,lat)
    #select a subset for checking the input
    dacheck=dacheck.sel(lon=valdata.lon,lat=valdata.lat)
    
    #we can correct the validation data for ignoring n=0, m=0 which has Y00=1 for all positions on the globe
    valdata=valdata-dain.basins*dain.time
    # # breakpoint()
    assert ((dacheck-valdata)/valdata).max().item() < tol


def test_synthesis_1d(generatedata):
    dain,valdata=generatedata
    #note: since valdata.lon and valdata.lat share the same dimensions analysis on the paris is performed rather than on a grid
    dacheck=dain.sh.synthesis(valdata.lon,valdata.lat)
    
    assert ((dacheck-valdata)/valdata).max().item() < tol

def test_synthesis_slice(generatedata,lon,lat):
    # test synthesis for the case where th input data is a view of the original memory
    dain,valdata=generatedata
    
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

