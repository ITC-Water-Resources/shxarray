# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#




from setuptools import setup,Extension
from setuptools_scm import get_version
from Cython.Build import cythonize
import Cython.Compiler.Options
from packaging.version import Version

import os 
import numpy as np
import sys

if sys.platform.startswith("win"):
    winplatform=True
    extra_args = ['/openmp']
else:
    winplatform=False
    extra_args = ['-fopenmp']

if "DEBUG_CYTHON" in os.environ:
    debug=True
    extra_args.append('-O0')
    #extra_args.append('-pg')
else:
    debug=False

#don't necessarily use cython
if "USE_CYTHON" in os.environ or winplatform:
    # note being on windows forces the use of cython
    useCython=True
    ext=".pyx"
    Cython.Compiler.Options.annotate = True
else:
    useCython=False
    ext=".cpp"

#Force the use of cython if numpy has a version < 2
if not useCython and Version(np.__version__) < Version ("2.0.0"):
    useCython=True
    ext=".pyx"


def listexts():
    names=["shlib"]
    exts=[]
    for nm in names:
        exts.append(Extension("shxarray."+nm.replace("/","."),["src/builtin_backend/"+nm+ext],include_dirs=[np.get_include(),"."], define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],extra_compile_args=extra_args,extra_link_args=extra_args))
    return exts

extensions=listexts()


if useCython:
    #additionally cythonize pyx files before building
    extensions=cythonize(extensions,language_level=3,annotate=True,gdb_debug=debug,compiler_directives={'embedsignature': True})

setup(
    version = get_version(root='.', relative_to=__file__),
    ext_modules=extensions
    )
