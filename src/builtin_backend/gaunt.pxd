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
    cdef cppclass Gaunt[T]:
        Gaunt() except +
        Gaunt(int n2,int n3,int m2, int m3) except +
        vector[T] get() except+
        T operator[](int j)except+
        int m()except+
        int nmin()except+
        int nmax()except+
    
    cdef cppclass GauntReal[T]:
        GauntReal() except +
        GauntReal(int n2,int n3,int mu2, int mu3) except +
        vector[T] get() except+
        int nmin()except+
        int nmax()except+
        int m()except+
        vector[pair[int,int]] nm() except+

# End of interface declaration



