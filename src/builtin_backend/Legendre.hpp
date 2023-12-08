/*! \file
 \brief 
 \copyright Roelof Rietbroek 2019
 \license
 This file is part of shxarray.
 */

#include<vector>
#include "cassert"
#ifndef FR_SH_LEGENDRE_HPP_
#define FR_SH_LEGENDRE_HPP_

#define UNNORM 0
#define GEODESY4PI 1


///@brief a class which computes and caches a (un)normalized Legendre polynomial
template<class ftype>
    class Legendre{
        public:
            Legendre(){}
            Legendre(int nmax):nmax_(nmax),pn_(nmax+1){}
            Legendre(int nmax,int norm):nmax_(nmax),pn_(nmax+1),norm_(norm){}
	    const std::vector<ftype> get(const ftype costheta);
        private:
            int nmax_=-1;
            std::vector<ftype> pn_{};
            int norm_=UNNORM;
    };


///@brief a class which computes and caches a (un)normalized Legendre polynomial
template<class ftype>
            const std::vector<ftype> Legendre<ftype>::get(const ftype costheta){
                assert(nmax_ >0);
                if (pn_[1] == costheta){
                    ///Quick return if already computed
                    return pn_;
                }
                
                ftype pnmin1=costheta;
                ftype pnmin2=1;
                ftype pn;
                pn_[0]=pnmin2;
                pn_[1]=pnmin1;

                for(int n=2;n<=nmax_;++n){
                   pn=((2*n-1)*costheta*pnmin1-(n-1)*pnmin2)/static_cast<ftype>(n);
                   pnmin2=pnmin1;
                   pnmin1=pn;
                   pn_[n]=pn;
                }
                
                /* possibly normalize */
                if (norm_ ==  GEODESY4PI){
                    for(int n=0;n<=nmax_;++n){
                        pn_[n]=pn_[n]/(2*n+1);

                    }
                }
                return pn_;
            }
            
#endif 
