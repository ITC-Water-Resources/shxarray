# Test spherical harmonic analysis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np
from shxarray.shlib import Ynm

nmax=20
tol=1e-7

#high degree testing
nmax=200
tol=1e-5

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



@pytest.fixture
def shinput():
    nbasins=5
    ntime=4
    basinscales=np.arange(nbasins)
    timescales=0.3*np.arange(ntime)

    scales1=xr.DataArray(basinscales,coords={"basins":basinscales})
    scales2=xr.DataArray(timescales,coords={"time":timescales})
    
    coords={"basins":basinscales,"time":timescales}
    dain_c=xr.DataArray.sh.ones(nmax ,auxcoords=coords,order='C')*scales1*scales2
    dain_f=xr.DataArray.sh.ones(nmax ,auxcoords=coords,order='F')*scales1*scales2
    return scales1,scales2,dain_c,dain_f

@pytest.fixture
def lon():
    return np.arange(-180,180.01,60.0)


@pytest.fixture
def lat():
    return np.arange(-60,60.01,60.0)

def test_analysis_full(shinput,validation,lon,lat):
    basinscales,timescales,dain_c,dain_f=shinput
    checkdata=validation
    #scale in the same way as the input data
    checkdata=checkdata*basinscales*timescales

    
    print("Checking C-storage order for full matrix")
    daout=dain_c.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    # breakpoint()
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

    print("Checking transposed C-storage order for full matrix")
    daout=dain_c.T.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol
    
    print("Checking F-storage order for full matrix")
    daout=dain_f.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

    print("Checking transposed F-storage order for full matrix")
    daout=dain_f.T.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

def test_analysis_slice(shinput,validation,lon,lat):
    basinscales,timescales,dain_c,dain_f=shinput
    dainsub_c=dain_c[1:4,2:4,1:]
    dainsub_f=dain_f[1:4,2:4,1:]
    
    checkdata=validation
    #note since we exclude SH coefficient n=0,m=0 the validation value can be corrected with 1.0
    checkdata=checkdata-1.0
    #scale in the same way as the input data
    checkdata=checkdata*basinscales[1:4]*timescales[2:4]

    print("Checking C-storage order for subset of matrix")
    daout=dainsub_c.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

    print("Checking transposed C-storage order for subset of matrix")
    daout=dainsub_c.T.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol
    
    print("Checking F-storage order for subset of matrix")
    daout=dainsub_f.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

    print("Checking transposed F-storage order for subset of matrix")
    daout=dainsub_f.T.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

def test_analysis_1d(shinput,validation,lon,lat):
    basinscales,timescales,dain_c,dain_f=shinput
    dain1d_c=dain_c[1,1,:]
    dain1d_f=dain_f[1,1,:]
    
    checkdata=validation
    #note since we exclude SH coefficient n=0,m=0 the validation value can be corrected with 1.0
    #scale in the same way as the input data
    checkdata=checkdata*basinscales[1].item()*timescales[1].item()

    print("Checking C-storage order for vector input")
    daout=dain1d_c.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol

    
    print("Checking F-storage order for subset of matrix")
    daout=dain1d_f.sh.analysis(lon,lat)
    dcheck=daout.sel(lon=checkdata.lon,lat=checkdata.lat)
    assert ((checkdata-dcheck)/dcheck).max().item() < tol


def test_lonlat_analysis(shinput,validation):

    basinscales,timescales,dain_c,dain_f=shinput
    checkdata=validation
    #scale in the same way as the input data
    checkdata=checkdata*basinscales*timescales

    print("Checking C-storage order for vector input")
    daout=dain_c.sh.analysis(checkdata.lon,checkdata.lat)
    # breakpoint()
    assert ((checkdata-daout)/daout).max().item() < tol

    
    print("Checking F-storage order for subset of matrix")
    daout=dain_f.sh.analysis(checkdata.lon,checkdata.lat)
    assert ((checkdata-daout)/daout).max().item() < tol
