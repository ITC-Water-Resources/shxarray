#  Xarray Extension for working with spherical harmonic data
This extension adds functionality to work with spherical harmonic data to xarray.


## Features and functionality (in progress)
* Gravity functionals: (convert from Stokes coefficients to various gravity functionals,s uch as equivalent water heights, geoid changes, etc.)
* Filters (e.g. Gaussian filters applied in the spectral domain)
* SNREI Earth models and functionals
* ...



## Installation
The extension makes use of cython generated code to speed up spherical harmonic synthesis and analysis.
**To install**:
1. (Optional, when a user installation is desired) create your own virtual environment with `venv` or Anaconda
2. clone this repository `git clone https://github.com/ITC-Water-Resources/shxarray.git`
3. Change to the repository directory  `cd shaxarray`
4. Install using pip  `pip install .`

**Install for development purposes**:
As above but use `pip install -e .` (editable install)
If cython code is being developed make sure to set the Enviroment variable:
`export USE_CYTHON=1`. This will regenerate the cython C/C++ code

Tip: when regenerating the shared library running `python setup.py build_ext` will be much faster than using pip. 


## Examples (in progress)
Jupyter notebook examples




