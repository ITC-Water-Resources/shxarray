# Test some basic operations of xarray objects filled with spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


import pytest
import shxarray
import xarray as xr
import os
import numpy as np
import time
from shxarray.core.logging import shxlogger
from shxarray.io.icgem import sample_tv_grav
import math

def test_icgem1_sig():
    """Test ICGEM Backends for xarray"""
    #note gzipped version (should also work)
    icgem1f=os.path.join(os.path.dirname(__file__),'testdata/icgem_test_sig_ITSG.gfc.gz')
    dsicgem=xr.open_dataset(icgem1f,engine='icgem')
    nmax=60
    assert dsicgem.sh.nmax == nmax, "nmax of icgem input is not 60"
    #note some data are missing so we only expect 599 entries
    assert dsicgem.sizes['nm'] == 599, "Not all dat has been read" 
    assert np.datetime64('2004-10-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcnm={(60,46):1.665788238155e-09}
    testsigcnm={(60,-46):2.363195963172e-12}
    for nm,cnm in testcnm.items():
        assert cnm == dsicgem.cnm.loc[:,nm].item(),f"Coefficient {nm} {cnm} not mathching"

    for nm,scnm in testsigcnm.items():
        assert scnm == dsicgem.sigcnm.loc[:,nm].item(),f"Sigma Coefficient {nm} {scnm} not mathching"


def test_icgem2_tv():
    """Test ICGEM2 Backends for xarray"""
    #note gzipped version (should also work)
    icgem1f=os.path.join(os.path.dirname(__file__),'testdata/CNES_GRGS.RL05MF_u202601.gfc.gz')
    nmax=91
    #also limit reading in higher than nmax 
    dsicgem=xr.open_dataset(icgem1f,engine='icgem',nmax=nmax)
    tsample=np.datetime64('2003-04-01') 
    dssampled=sample_tv_grav(dsicgem,tsample)

    assert dssampled.sh.nmax == nmax, f"nmax of icgem input is not {nmax}"

   # gfct  90   87  1.90973490376E-09  3.06261539611E-10 7.2870E-15 7.8750E-15 20020416.1016 20040101.0000
# trnd  90   87 -4.80147817178E-14 -4.57662031755E-14 1.2430E-14 1.3550E-14 20020416.1016 20040101.0000
# acos  90   87  2.26861718802E-14  1.62981340831E-14 3.7150E-15 3.9460E-15 20020416.1016 20040101.0000 1.0
# asin  90   87 -2.46353148623E-15 -5.52179488317E-15 2.8390E-15 3.1620E-15 20020416.1016 20040101.0000 1.0
# acos  90   87 -7.27574781068E-15 -2.60867498328E-15 2.9610E-15 3.2430E-15 20020416.1016 20040101.0000 0.5
# asin  90   87  1.39373205841E-15 -5.34176729895E-15 3.3680E-15 3.5410E-15 20020416.1016 20040101.0000 0.5 

# gfc   91   89  2.65649920873E-09 -3.52571733178E-09 3.8550E-12 3.8560E-12
    t0=np.datetime64('2002-01-01')
    # t0=np.datetime64('2002-04-16T10:16')
    dt_yr=(tsample-t0)/np.timedelta64(int(365.25*24),'h')
    ohm=2*np.pi*dt_yr
    testcnm={(90,87):1.90973490376E-09+(-4.80147817178E-14*dt_yr)+np.cos(ohm)* 2.26861718802E-14 + np.sin(ohm)*-2.46353148623E-15 + np.cos(2*ohm)*(-7.27574781068E-15 )+np.sin(2*ohm)*(1.39373205841E-15 )}
    testcnm[(0,0)]=1
# gfc  300  300 -1.25298319912E-11 -1.00609822622E-11 4.0920E-11 4.0920E-11
    testcnm[(91,-89)]=-3.52571733178E-09 
    
    for nm,cnm in testcnm.items():
        assert math.isclose(cnm,dssampled.cnm.loc[:,nm].item()),f"Coefficient {nm} {cnm} not matching"



def test_icgem_nosig():
    """Test ICGEM Backend for xarray for file without sigmas"""
    icgem2f=os.path.join(os.path.dirname(__file__),'testdata/icgem_test_nosig_ITSG.gfc')
    dsicgem=xr.open_dataset(icgem2f,engine='icgem')
    nmax=2
    assert dsicgem.sh.nmax == nmax, "nmax of icgem input is not 2"
    assert dsicgem.sizes['nm'] == 9, "Not all dat has been read" 
    assert np.datetime64('2013-11-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcnm={(2,0):-4.841696152100e-04,(1,1):0.0}
    for nm,cnm in testcnm.items():
        assert cnm == dsicgem.cnm.loc[:,nm].item(),f"Coefficient {nm} {cnm} not mathching"


def test_icgem_static_tv():
    """Test ICGEM Backend for xarray for file with static and parametric coefficients"""
    icgem3f=os.path.join(os.path.dirname(__file__),'testdata/ITSG-Grace2018s.gfc.gz')
    dsicgem=xr.open_dataset(icgem3f,engine='icgem')
    nmax=200
    assert dsicgem.n.max().item(), "nmax of icgem input is not 200"
    # assert dsicgem.sizes['nm'] == 9, "Not all dat has been read" 
    # assert np.datetime64('2013-11-15T12:00:00.000000000') ==dsicgem.time.data[0], "Time epoch is off or not present"
    testcoef={('gfct',(10,2)):-9.396504743590e-08,('gfct',(200,-200)):1.956728673164e-10,
              ('trnd',(195,159)):7.843896315171e-17,('trnd',(195,-159)):4.073388066625e-17,
              ('acos', (7,3)): 2.269912528926e-11,('acos', (7,-3)): -3.872126870042e-12} 
    
    testsigcoef={('trnd',(195,159)):7.023988334695e-15,('trnd',(195,-159)):7.022127457848e-15, ('acos',(7,-3)): 4.115793871821e-13,('acos',(7,3)): 4.108227520945e-13}

    # acos    7    3  2.269912528926e-11 -3.872126870042e-12  4.108227520945e-13  4.115793871821e-13 1.0
    for (var,nm),cnm in testcoef.items():
        assert cnm == dsicgem.cnm.sel(ityp=var).sel(n=nm[0]).sel(m=nm[1]).item(),f"Coefficient {nm} {cnm} not matching"

    for (var,nm),cnm in testsigcoef.items():
        assert cnm == dsicgem.sigcnm.sel(ityp=var).sel(n=nm[0]).sel(m=nm[1]).item(),f"Sigma Coefficient {nm} {cnm} not matching"


def test_gsm_cnes():
    """Test GSM Backends (CNES version) for xarray"""
    #note gzipped version (should also work)
    cnesf=os.path.join(os.path.dirname(__file__),'testdata/GSM-2_2023121-2023151_GRFO_CNESG_TSVD_0500.txt.gz')
    dscnes=xr.open_dataset(cnesf,engine='gsmv6')
    nmax=90
    assert dscnes.sh.nmax == nmax, f"nmax of icgem input is not {nmax}"
    #note some data are missing so we only expect 599 entries
    assert dscnes.sizes['nm'] == (nmax+1)*(nmax+2)-nmax-1, "Not all dat has been read" 
    assert np.datetime64('2023-05-16T00:00:00.000000000') ==dscnes.time.data[0], "Time epoch is off or not present"
    testcnm={(90,90):0.881655402494E-09}
    testsigcnm={(90,-90): 9.6427E-15}
    
    testcnm[(90,-83)]=-0.120965152481E-09
    testsigcnm[(90,83)]: 1.6863E-15

    for nm,cnm in testcnm.items():
        assert cnm == dscnes.cnm.loc[:,nm].item(),f"Coefficient {nm} {cnm} not mathching"

    for nm,scnm in testsigcnm.items():
        assert scnm == dscnes.sigcnm.loc[:,nm].item(),f"Sigma Coefficient {nm} {scnm} not mathching"
