# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++


cimport cython
from libcpp.vector cimport vector
from libcpp.pair cimport pair
# C++ / Cython interface declaration 


cdef extern from "Gaunt.hpp":
    cdef cppclass Gaunt[T] nogil:
        Gaunt() except +
        Gaunt(int n2,int n3,int m2, int m3) except +
        void set (int n2,int n3,int m2, int m3) nogil
        vector[T] get() except+
        T operator[](int j)except+
        int m()except+
        int nmin()except+
        int nmax()except+
    
    cdef cppclass GauntReal[T] nogil:
        GauntReal() except +
        GauntReal(int nmax) nogil
        GauntReal(int n2,int n3,int mu2, int mu3) nogil
        void set(int n2,int n3,int mu2, int mu3) nogil
        T operator[](cython.size_t i)
        T at(int n1,int m1) except+
        vector[T] get() except+
        cython.size_t size() 
        int nmin()except+
        int nmax()except+
        vector[pair[int,int]] & nmvec() except+
        pair[int,int] & nm(cython.size_t i)
# End of interface declaration



