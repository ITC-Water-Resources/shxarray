# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#




from setuptools import setup,Extension
from Cython.Build import cythonize
import Cython.Compiler.Options
import os 
import numpy as np

debug=False
#don't necessarily use cython
if "USE_CYTHON" in os.environ:
    useCython=True
    ext=".pyx"
    Cython.Compiler.Options.annotate = True
else:
    useCython=False
    ext=".cpp"


def listexts():
    names=["shlib"]
    exts=[]
    for nm in names:
        exts.append(Extension("shxarray."+nm.replace("/","."),["src/builtin_backend/"+nm+ext],include_dirs=[np.get_include(),"."], define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]))
    return exts

extensions=listexts()


if useCython:
    #additionally cythonize pyx files before building
    extensions=cythonize(extensions,language_level=3,annotate=True,gdb_debug=debug)

setup(
    ext_modules=extensions
    # entry_points={
        # "xarray.backends": [
            # "icgem=shxarray.io.shiobackend:ICGEMBackEntryPoint",
            # "gsmv6=shxarray.io.shiobackend:GSMv6BackEntryPoint"
        # ],
        # "shxarray.computebackends": [
            # "builtin=shxarray.shlib:SHComputeBackend",
            # ],
    # }
    )
