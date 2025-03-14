from cython cimport size_t
from libcpp.pair cimport pair
cdef extern from "Helpers.hpp":
    cdef cppclass Nmindex nogil:
        Nmindex() except +
        Nmindex(size_t maxsize) except +
        void set(pair[int,int] nm, size_t ix)
        size_t operator[] (const pair[int,int] & nm)
