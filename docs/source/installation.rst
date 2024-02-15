Installation and General Usage
==============================
The latest **shxarray** is hosted on `pypi <https://pypi.org/project/shxarray/>`_ and can be installed through pip: 

``pip install shxarray`` 

Part of the module is written is `Cython <https://cython.readthedocs.io/en/latest/>`_, which means that a c compiler is needed to build the package. A binary wheel is currently not offered, but this may be offered in the future.

Import and usage
----------------
For most operations, a simple import will expose the xarray extensions. For example:

.. code-block:: python
   
   import shxarray
   import xarray as xr

   #Initialize a dataarray with zeros which has a dimension spanning degrees from nmin to nmax

   nmax=20
   nmin=2
   dazeros=xr.DataArray.sh.ones(nmax=nmax,nmin=nmin)



Development
-----------
Users interested in developing can install the latest version from `github <https://github.com/ITC-Water-Resources/shxarray/tree/main>`_. Cython is needed in case the binary extension is being developed, and users can consult the dedicated instructions on the github repository.

Code can be supplied with `numpy docstrings <https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html>`_ so they can be parsed into this documentation.


