# Test some basic operations of xarray objects filled with spherical harmonic datasets
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023


import pytest
import shxarray
import xarray as xr
import os
import numpy as np

def test_icgem():
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


@pytest.fixture
def sinexval():
    # some validation values based on prior visual inspection of sinex file
    #+SOLUTION/ESTIMATE
    # ...
    # 91 CN        9 --    7 03:259:00000 ---- 2 -1.17979570334447e-07 1.01491e-12
    # 92 SN        9 --    7 03:259:00000 ---- 2 -9.69270369824925e-08 1.02558e-12
    # ..
    #+SOLUTION/APRIORI
    #...
    # 25 CN        5 --    2 03:259:00000 ---- 2  6.52120740523400e-07 0.00000e+00
    # 26 SN        5 --    2 03:259:00000 ---- 2 -3.23349434999185e-07 0.00000e+00
    #..
    #+SOLUTION/NORMAL_EQUATION_VECTOR
    #*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___RIGHT_HAND_SIDE___
    #1 CN        2 --    0 03:259:00000 ---- 2 -3.27811136401904e+12
    #2 CN        2 --    1 03:259:00000 ---- 2 -4.32367855180938e+11
    #+SOLUTION/NORMAL_EQUATION_MATRIX U
    # ....
    # 7     7  2.92598076849103e+23 -3.55351048890917e+21 -3.78432268065340e+21
    # 7    10 -6.34992891236174e+21  2.34638506385098e+20  4.78761459861900e+21
    # ...
    valdata={} 
    valdata["sol_est"]=[(9,7, -1.17979570334447e-07),(9,-7,-9.69270369824925e-08)] 
    valdata["sol_std"]=[(9,7,1.01491e-12),(9,-7, 1.02558e-12)]
    valdata["apri_est"]=[(5,2,6.52120740523400e-07),(5,-2,-3.23349434999185e-07)] 
    valdata["rhs"]=[(2,0,-3.27811136401904e+12),(2,1,-4.32367855180938e+11)]
    valdata["N"]=[(7,7, 2.92598076849103e+23),(7,8,-3.55351048890917e+21),(7,9,-3.78432268065340e+21),(7,10,-6.34992891236174e+21),(7,11,2.34638506385098e+20),(7,12,4.78761459861900e+21)]
    return valdata

     


def test_sinex(sinexval):
    """ Test reading of Normal equation systems in SINEX format"""
    import requests
    url="https://ftp.tugraz.at/outgoing/ITSG/GRACE/ITSG-Grace2018/monthly/normals_SINEX/monthly_n96/ITSG-Grace2018_n96_2003-09.snx.gz"
    sinexfile=os.path.join(os.path.dirname(__file__),'testdata', os.path.basename(url))
    
    if not os.path.exists(sinexfile):
        print(f"Downloading {sinexfile}...")
        r=requests.get(url)
        with open(sinexfile,'wb') as fid:
            fid.write(r.content)
    # sinexfile=os.path.join(os.path.dirname(__file__),'testdata/ITSG-Grace2018_n96_2003-09.snx') 
    #quick read which stops when encountering a normal matrix
    #Engine does not need to be specified because file corresponds to commonly used filename pattern for sinex 
    dsneqsinex=xr.open_dataset(sinexfile,drop_variables=["N"])

    for var in ["sol_est","apri_est","sol_std","rhs"]:
        for n,m,val in sinexval[var]:
            assert val == dsneqsinex[var].sel(n=n,m=m).item()
        # assert dsneqsinex[var].sel(n=
    
    #read version with entire normal matrix
    dsneqsinex=xr.open_dataset(sinexfile,engine='sinex')
    for ix,iy,val in sinexval["N"]:
        #note index ix,and iy are 1-indexed
        assert val == dsneqsinex.N[ix-1,iy-1].item()



