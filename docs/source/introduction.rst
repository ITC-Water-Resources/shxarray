Introduction
============

The **shxarray** package is aimed to make using spherical harmonic operations more accessible to a community which is used to xarray.

Putting degrees and orders in a MultiIndex
------------------------------------------

Spherical harmonic coefficients are often put in 2-dimensional arrays with degree and order spanning the 2 dimensions. This has the advantage that individual coefficients can be easily referenced as e.g. `cnm=mat[n,m]`. However, since only the upper triangle of those matrices are non-zero, the sparseness of these matrices is not made use of. When working with large datasets which also span other dimensions such as time or different levels, this will cause large segments of zeros.

A `pandas.MultiIndex <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.MultiIndex.html>`_ can facilitate the stacking of degrees and orders in a single dimension. On top of that, this multindex can then be used in `xarray` to work with spherical harmonics along a single coordinate. In **shxarray**, the spherical harmonic index is generally denoted as ``nm``, while when two spherical harmonic coordinates are needed alternative versions such as ``nm_`` are added. 

Exposing spherical harmonic functionality to xarray
---------------------------------------------------

Some of the functionality needed for working with spherical harmonics need specialized access to the degree and order information. The aim of **shxarray** is to expose common functionality through `xarray accessors <https://docs.xarray.dev/en/stable/internals/extending-xarray.html>`_. This allows, besides familiar syntax for xarray users, also chaining of operations in a compact syntax.

Delegate common operations to xarray
------------------------------------

In contrast to specialized spherical harmonic operations, many operations on spherical harmonic data can be delegated to xarray itself. Wherever possible, functionality and broadcasting features of xarray are made use for a consistent syntax.





