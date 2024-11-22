# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#
# distutils: language = c++


cimport cython
from libcpp.vector cimport vector
# C++ / Cython interface declaration 
cdef extern from "Wigner3j.hpp":
    cdef cppclass Wigner3j[T] nogil:
        Wigner3j() except +
        Wigner3j(int j2,int j3,int m2, int m3) except +
        vector[T] get() except+
        void set(int j2,int j3,int m2, int m3) except+
        T operator[](int j)except+
        int jmin()
        int jmax()
        int m()

# End of interface declaration



