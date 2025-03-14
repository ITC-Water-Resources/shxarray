import xarray as xr
import shxarray
from math import sqrt
import numpy as np

def p2s(daobj):
    #Initialize the output array
    nmax=int(daobj.sh.nmax/2)
    dsp2s=xr.DataArray.sh.zeros(nmax,nshdims=2)
    nsh=len(dsp2s['nm'])
    ntot=0
    norm=sqrt(4*np.pi)
    for i in range(nsh):
        n2,m2=dsp2s.nm[i].item() 
        for j in range(i+1):
            n3,m3=dsp2s.nm_[j].item() 
            dagaunt=xr.DataArray.sh.gauntReal(n3,n2,m3,m2)
            ntot+=1
            dsp2s.data[i,j]=norm*daobj.dot(dagaunt).item()
            if i != j:
                #also mirror data
                dsp2s.data[j,i]=dsp2s.data[i,j]
    return dsp2s
