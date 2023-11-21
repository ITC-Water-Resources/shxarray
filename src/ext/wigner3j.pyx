# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++
# cython: profile=False

cimport cython
from wigner3j cimport Wigner3j

from cython.operator cimport dereference as deref
# import numpy as np 
cimport numpy as np
from libc.math cimport sqrt

# Cython wrapper class
cdef class W3j:
    """Double precision Wigner 3J symbols"""
    cdef Wigner3j[double]*w3j_ptr  # Pointer to the wrapped C++ class
    def __cinit__(self, int j2,int j3,int m2, int m3):
        self.w3j_ptr = new Wigner3j[double](j2,j3,m2,m3)
    
    def __dealloc__(self):
        del self.w3j_ptr

    def __getitem__(self,int j):
        return deref(self.w3j_ptr)[j]

    property data:
        def __get__(self):
            return np.asarray(deref(self.w3j_ptr).get())

