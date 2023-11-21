/*! \file
 \brief 
 \copyright Roelof Rietbroek 2023
 \license
 This file is part of shxarray.
 */


#include <cmath>
#ifndef HELPERS_HPP_
#define HELPERS_HPP_


inline int csphase(int m){
    return 1-2*(std::abs(m)%2);
}

inline int kronecker(int n1,int n2){
   return n1==n2?1:0;
}

#endif 
