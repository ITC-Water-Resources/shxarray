# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++


cimport cython
from libcpp.vector cimport vector
from libcpp.pair cimport pair
# C++ / Cython interface declaration 
cdef extern from "Legendre.hpp":
    cdef cppclass Legendre[T] nogil:
        Legendre(int nmax) except +
        Legendre(int nmax,int norm) except +
        vector[T] get(T costheta) except+

# C++ associated Legendre functions
cdef extern from "Legendre_nm.hpp":
    cdef cppclass Legendre_nm[T] nogil:
        Legendre_nm() except +
        Legendre_nm(int nmax) except +
        void set(T costheta, double arr[] ) except+
        @staticmethod
        cython.size_t i_from_nm(int n,int m, int nmax)
        cython.size_t idx(int n,int m)
        @staticmethod
        pair[int,int] nm_from_i(cython.size_t idx, int nmax)
        #cached version:
        pair[int,int] nm(cython.size_t idx)
        int nmax()
        cython.size_t size()
# End of interface declaration

