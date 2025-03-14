# Test Wigner 3j symbols and Gaunt  
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023



import pytest
import shxarray
import xarray as xr
import numpy as np
import gzip
import os


@pytest.fixture
def product2sumVal():
    # """Create a on-disk dataset with product 2 sum matrix based on ocean coefficients"""
    nmax=10
    p2sumfile=os.path.join(os.path.dirname(__file__),f'testdata/P2Sum_ocean{nmax}.nc')
    if not os.path.exists(p2sumfile):
        print(f"Generating {p2sumfile} using rlftlbx")
        import requests
        import subprocess
        from shxarray.io.binv_legacy import readBINV
        oceurl="https://github.com/strawpants/geoshapes/raw/refs/heads/master/raster/ocean/ne_10m_oceansh_n300.nc"
        ocet=os.path.join("/tmp",os.path.basename(oceurl))
        if not os.path.exists(ocet):
            with requests.Session() as s:
                r=s.get(oceurl)
                with open(ocet,'wb') as fid:
                    fid.write(r.content)
        from io import StringIO,BytesIO  
        #Load data and rewrite to a file with 
        dsoce=xr.open_dataset(ocet).sh.truncate(2*nmax)
        ocesht='/tmp/oce.sh'
        with open(ocesht,'wt') as fid:
            oceascii=dsoce.oceansh.sh.to_ascii(fid)

        #Create the validation file using the RLFtlbx tool box
        p2sumtmp='/tmp/p2sum.bin'
        with open(p2sumtmp,'wb') as fid:
            outp=subprocess.run(['SH_prod2sum_mat_openmp',f"-l{nmax}",ocesht],capture_output=True)
            fid.write(outp.stdout)

        dsp2s=readBINV(p2sumtmp)

        #add original ocean function coefficients
        dsp2s['oceansh']=dsoce.rename(nm='nm_orig',n='n_orig',m='m_orig').oceansh
        #save to the p2sumfile
        dsp2s.reset_index(['nm','nm_']).to_netcdf(p2sumfile)
    else:
        dsp2s=xr.open_dataset(p2sumfile)


    return dsp2s



def test_product2sum(product2sumVal):
    rtol=1e-8
    #generate Product to sum matrix (note result will have a nmax which is half the size of the input max degree
    dsocean=product2sumVal.oceansh.rename(nm_orig='nm',n_orig='n',m_orig='m').sh.build_nmindex()
    
    dsp2s=dsocean.sh.p2s().T
    daval=product2sumVal.mat.sh.build_nmindex().sh.build_nmindex('_').loc[dsp2s.nm,dsp2s.nm_] 
    closeEnough=np.allclose(daval.data,dsp2s.data,rtol=rtol)
    
    # import pdb;pdb.set_trace()
    assert(closeEnough)    




