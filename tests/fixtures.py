from xarray.core.dataarray import DataArray
from shxarray.kernels import ParabolicCap, Disk
import pytest
import xarray as xr
import os
import numpy as np
import shxarray

nmax=200

@pytest.fixture
def shcapsvalidation():
    nmax=200
    lon=[135 , -71]
    lat=[-31, 65]
    radius=10 #radius of parabolic cap in deg
    pcap=ParabolicCap(nmax,radius)
    pcapsh=pcap.position(lon,lat)
    return pcapsh.rename(dict(lon="plon",lat="plat",nlonlat="npoints"))


@pytest.fixture
def shcapsvalidationaux(shcapsvalidation):
    # add another set of dimension

    nbasins=5
    ntime=4
    basinscales=np.arange(nbasins)
    timescales=0.3*np.arange(ntime)

    shcapaux=shcapsvalidation.expand_dims({"time":len(timescales),"basins":len(basinscales)}).assign_coords(dict(time=timescales,basins=basinscales))
    shcapaux=shcapaux*shcapaux.basins*shcapaux.time   
    return shcapaux


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
    return dagrd.rename(time="npoints")



# The validation datatset is produced from the RLFTLBX using the following commands:
# SH_dump -l200 -u | SH_2_nc  -R=-180/180/-60/60 -I=60 -g
# Which creates a unit spherical harmonic dataset mapped to a global grid
@pytest.fixture
def grdvalidation():
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
def generategrddata(grdvalidation,request):
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
    valdata=grdvalidation*scales1*scales2

    return dain,valdata

@pytest.fixture
def basin_sim_data():
    """
    Generate an simulated (annual/semiannual) signal in three non-overlapping circular basins, where amplitude and phase and size of the basins differ
    """
    time= np.arange(np.datetime64("2000-01-01"), np.datetime64("2010-12-31"), np.timedelta64(7, "D"))
    
    #generate a set of circular basins which are closeby but are non overlapping
    nmax=120
    radius=[11.5,9.5,13.5] #radii in degrees of the disks
    lon=[-51.38,-63,-77.34]
    lat=[-10.22,8,-10.30]
    dabasins=None
    for rad,lo,la in zip(radius,lon,lat):
        pdisk=Disk(nmax*2,rad)
        
        datmp=pdisk.position([lo],[la]).rename(dict(lon="plon",lat="plat",nlonlat="basins"))
        if dabasins is None:
            dabasins=datmp
        else:
            dabasins=xr.concat([dabasins,datmp],dim="basins")
    # create a seasonal signal in the basins and a phase shifted seasonal signal with lower amplitude in the surrounding areas 
    # amplitudes 
    amplitudes=np.array([3.,20.,1.])
    phases=np.array([0,90/365,200/365])

    semiamplitudes=0.1*amplitudes
    semiphases=phases+1/12

    #generate a seasonal signal
    t0=np.datetime64("2005-01-01")
    time_yrs=((time-t0)/np.timedelta64(365,'D')).astype(float)
    seas_cos=np.cos(2*np.pi*(np.matlib.repmat(time_yrs,len(phases),1).T+np.matlib.repmat(phases,len(time_yrs),1)))
    seas_cos*=np.matlib.repmat(amplitudes,len(time_yrs),1)

    #generate a semiannual signal
    semiseas_cos=np.cos(4*np.pi*(np.matlib.repmat(time_yrs,len(semiphases),1).T+np.matlib.repmat(semiphases,len(time_yrs),1)))
    semiseas_cos*=np.matlib.repmat(semiamplitudes,len(time_yrs),1)

    #add the contribution of the basins together
    datws=(xr.DataArray(seas_cos+semiseas_cos,dims=['time','basins'],coords=dict(time=time))*dabasins).sum('basins')
    
    ds=datws.to_dataset(name="tws")
    ds['basin_sh']=dabasins
    #add truth-mean series
    ds['basin_avs']=(["time","basins"],seas_cos+semiseas_cos)
    return ds
     




