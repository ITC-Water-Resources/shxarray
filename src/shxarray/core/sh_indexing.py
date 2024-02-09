# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import pandas as pd
from enum import IntEnum
from functools import total_ordering


@total_ordering
class trig(IntEnum):
    c=0
    s=1


class SHindexBase:
    # _nmax=None
    # _nmin=None
    # _squeeze=True
    # def __init__(self,nmax=None,nmin=None,squeeze=True):
        # _nmax=nmax
        # _nmin=nmin
        # _squeeze=squeeze

    @staticmethod
    def nsh(nmax,nmin=0,squeeze=True):
        """
        Compute the total amount of spherical harmonic coefficients
        :param nmax: maximum spherical harmonic degree
        :param nmin: minimum spherical harmonic degree
        :param squeeze: When true, don't store zero order Sine coefficients (will be zero)
        :return Amount of spherical harmonic coefficients
        """
        assert nmax>=nmin

        sz=(nmax+1)*(nmax+1)
        if not squeeze:
            sz+=nmax+1

        if nmin > 0:
            #possibly remove the number of coefficients which have n < nmin (calls this function itself)
            sz-=SHindexBase.nsh(nmin-1,0,squeeze=squeeze)

        return sz

    @staticmethod
    def nmt_mi(nmax,nmin=0,squeeze=True):
        """ create a multindex guide which varies with n, then m, and than trigonometric sign"""
        if squeeze:
            nmt=[(n,m,t) for t in [trig.c,trig.s] for n in range(nmin,nmax+1) for m in range(n+1) if not (m == 0 and t == trig.s) ]
        else:
            nmt=[(n,m,t) for t in [trig.c,trig.s] for n in range(nmin,nmax+1) for m in range(n+1)]
        return SHindexBase.mi_fromtuples(nmt)

    @staticmethod
    def shi(nmax,nmin=0,squeeze=True,dim="shi"):
        """Convenience function which returns a dictionary which can be used as input for xarray constructors"""
        return {dim:(dim,SHindexBase.nmt_mi(nmax,nmin,squeeze))}


    @staticmethod
    def mi_fromtuples(nmt):
        return pd.MultiIndex.from_tuples(nmt,names=["n","m","t"])
    
    @staticmethod
    def mi_fromarrays(nmt):
        return pd.MultiIndex.from_arrays(nmt,names=["n","m","t"])
    
    @staticmethod
    def mi_toggle(mi):
        """Rename the levels of the multindex so that they can be use as transposed versions"""
        if "n" in mi.names:
            return mi.rename([nm+"_" for nm in mi.names])
        else:
            return mi.rename([nm[:-1] for nm in mi.names])

