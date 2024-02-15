# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import pandas as pd


class SHindexBase:
    name="nm"
    name_t="nm_"
    @staticmethod
    def nsh(nmax,nmin=0, squeeze=True):
        """
        Compute the total amount of spherical harmonic coefficients for a given range

        Parameters
        ----------
        nmax : int
            maximum spherical harmonic degree
        nmin : int, optional
            minimum spherical harmonic degree
        squeeze: bool,optional
            Legacy option used when Sine coefficients which have m=0 need to be included
            

        Returns
        -------
        int
            The amount of spherical harmonic coefficients in this range
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
    def nm_mi(nmax,nmin=0):
        """
        Generate a MultiIndex of degree and order which span a spherical harmonic degree range
        
        In the case of real spherical harmonic, orders m < 0  denote Sine coefficients

        Parameters
        ----------
        nmax : int
            maximum spherical harmonic degree
        nmin : int, optional
            minimum spherical harmonic degree

        Returns
        -------
        pandas.MultiIndex
            A MultiIndex with degrees "n" and orders "m" 
        """
        nm=[(n,m) for n in range(nmin,nmax+1) for m in range(-n,n+1)]
        return SHindexBase.mi_fromtuples(nm)

    @staticmethod
    def nm(nmax,nmin=0):
        """
        Convenience function which returns a dictionary which can be used as input for xarray constructors
        
        Parameters
        ----------
        nmax : int
            maximum spherical harmonic degree
        nmin : int, optional
            minimum spherical harmonic degree

        Returns
        -------
        dictionary
            A dictionary specifying the degree and orders and corresponding dimension names
            in the form of {dim:(dim,nm)}
        """
        return {SHindexBase.name:(SHindexBase.name,SHindexBase.nm_mi(nmax,nmin))}


    @staticmethod
    def mi_fromtuples(nm):
        """
        Generate a MultiIndex of degree and order from a list of (degree,order) tuples
        
        In the case of real spherical harmonic, orders m < 0  denote Sine coefficients

        Parameters
        ----------
        nm : list
            A list of tuples with degree and order

        Returns
        -------
        pandas.MultiIndex
            A MultiIndex with degrees "n" and orders "m" 
        """

        return pd.MultiIndex.from_tuples(nm,names=["n","m"])
    
    @staticmethod
    def mi_fromarrays(nm):
        """
        Generate a MultiIndex of degree and order from an array of degree and order [[n..],[..m]]
        
        In the case of real spherical harmonic, orders m < 0  denote Sine coefficients

        Parameters
        ----------
        nm : array-like
            An array which hold a vector of degrees and orders

        Returns
        -------
        pandas.MultiIndex
            A MultiIndex with degrees "n" and orders "m" 
        """
        return pd.MultiIndex.from_arrays(nm,names=["n","m"])
    
    @staticmethod
    def mi_toggle(mi,ending=''):
        """
        Rename the levels of a (nm)-multindex so that they can be use as alternative coordinates (e.g. transposed versions)

        The levels will be swicthed back and fort between the following formats
        oldname <-> oldname_[ending]

        Parameters
        ----------
        mi : pandas.MultiIndex
            A MultiIndex with degree and orders
        ending: str, optional
            A string which can be additionally appended

        Returns
        -------
        pandas.MultiIndex
            A MultiIndex with renamed levels 
        """

        app="_"+ending
        applen=len(app)
        if "n" in mi.names:
            return mi.rename([nm+app for nm in mi.names])
        else:
            return mi.rename([nm[:-applen] for nm in mi.names])

