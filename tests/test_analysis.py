# Test spherical harmonic analysis
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025



from fixtures import shcapsvalidation as validation
from fixtures import grdinput


tol=1e-7

def test_analysis_shlib(grdinput,validation):
    """ Test the SH analysis to see if an acceptable sh is computued"""
    nmax=validation.sh.nmax
    dsgrd=grdinput
    checkdata=validation
    dscheck=dsgrd.sh.analysis(nmax)
    
    #ok we need to rename the input dimension time
    dadiff=checkdata-dscheck#.rename(time='npoints')
    print("Checking whether retrieved SH solution is within tolerance")
    maxdiff=max(abs(dadiff.max()),abs(dadiff.min())) 
    assert(maxdiff < tol)



