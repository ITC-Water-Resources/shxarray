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

from shxarray.io.sinex import read_sinex

def test_sinex():
    import requests
    url="https://ftp.tugraz.at/outgoing/ITSG/GRACE/ITSG-Grace2018/monthly/normals_SINEX/monthly_n96/ITSG-Grace2018_n96_2003-09.snx.gz"
    sinexfile=os.path.join(os.path.dirname(__file__),'testdata', os.path.basename(url))
    
    if not os.path.exists(sinexfile):
        print(f"Downloading {sinexfile}")
        r=requests.get(url)
        with open(sinexfile,'wb') as fid:
            fid.write(r.content)

    dsneqsinex=read_sinex(sinexfile)    



