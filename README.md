#  Xarray Extension for working with spherical harmonic data
[![Build](https://github.com/ITC-Water-Resources/shxarray/actions/workflows/python-publish.yml/badge.svg)](https://github.com/ITC-Water-Resources/shxarray/actions/workflows/python-publish.yml)
[![PyPI version](https://badge.fury.io/py/shxarray.svg)](https://badge.fury.io/py/shxarray)
[![Documentation Status](https://readthedocs.org/projects/shxarray/badge/?version=latest)](https://shxarray.wobbly.earth/latest/?badge=latest)

This extension adds functionality to work with Spherical Harmonics to [Xarray](https://github.com/pydata/xarray).

## Features and functionality (in progress)
* Gravity functionals: (convert from Stokes coefficients to various gravity functionals, such as equivalent water heights, geoid changes, etc.)
* Filter (e.g. Gaussian or anisotropic filter applied in the spectral domain)
* The use of Xarray like operations allow for applying functionality to multi-dimensional datasets


## Getting started
The tutorials in the [documentation](https://shxarray.wobbly.earth/stable/tutorial.html) will gradually epxand and provide Jupyter notebook workflows to start with.


## Installation
You can install this package from PyPi using:
```
pip install shxarray
```


## Development Installation
If you want to help in the development of this package, it's best to clone the repository to allow for modifications and pull requests. The extension makes use of [Cython](https://cython.readthedocs.io/en/latest/) generated code to speed up spherical harmonic synthesis and analysis.

1. Create your own virtual environment with `venv` or Anaconda *(Optional but recommended, when a user installation is desired)*
2. Clone this repository `git clone https://github.com/ITC-Water-Resources/shxarray.git`
3. Change to the repository directory `cd shaxarray`
4. Set the environment variable `export USE_CYTHON=1` *(Optional and only in the case Cython code is being developed or needs to be regenerated)*
5. Install using pip  `pip install .` or use `pip install -e .` for an editable install
 
### Cython build tip on an editable install
From the repository root directory, regenerating the shared library running 

```python ./setup.py build_ext``` 

will be much faster than using 

```pip install -e .``` 


This will build the shared library in for example `./build/lib.linux-x86_64-cpython-3xx/shxarray/shlib.cpython-3xx-x86_64-linux-gnu.so`. To make sure changes are picked up in your editable install you can create a symbolic link in the Python part of the library e.g. :

```
cd src/shxarray/
ln -sf ../../build/lib.linux-x86_64-cpython-311/shxarray/shlib.cpython-311-x86_64-linux-gnu.so
```


## Contributing
This repository is under development and contributions and feedback are welcome.

### Contributors
* Main developer: Roelof Rietbroek (r.rietbroek@utwente.nl)
* ..



