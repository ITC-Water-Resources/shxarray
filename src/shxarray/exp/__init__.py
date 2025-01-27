"""
Experimental compute backend for shxarray

This submodule can be used to test and stage functionality in pure Python before moving it to a `shlib compute engine <shxarray.shlib.html>`_.

Functionality in this submodule is available when using the ``engine='exp'`` option in the xarray accessor routines.

Developer Notes:
    To add functionality to the engine add a Python file in the ``exp`` directory and import it in the ``__init__.py`` file.

"""
from shxarray.exp.p2s import *
