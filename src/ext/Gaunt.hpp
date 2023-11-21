/*! \file
 \brief 
 \copyright Roelof Rietbroek 2023
 \license
 This file is part of shxarray.
 */

#include<vector>
#include "cassert"
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <iostream>
#include "Helpers.hpp"

#ifndef FR_SH_GAUNT_HPP_
#define FR_SH_GAUNT_HPP_





///@brief class to compute Gaunt coefficients (integral of 3 (conplex) spherical harmonic base functions)
template<class ftype>
class Gaunt{
    public:
        Gaunt(){};
        Gaunt(int n2,int n3,int m2,int m3);
        std::vector<ftype> get ()const{return data;}
        ftype operator[](int n)const{return data[idx_.at(n)];}
        int nmin()const{return n1min_;}
        int nmax()const{return n1max_;}
        int m()const{return m1_;}
    private:
        int n1max_=-1;
        int n1min_=-1;
        int m1_=-1;
        int n2_=-1;
        int m2_=-1;
        int n3_=-1;
        int m3_=-1;
        size_t sz_=0;
        std::vector<ftype> data;
        std::map<int,size_t> idx_;
};



template<class ftype>
Gaunt<ftype>::Gaunt(int n2,int n3, int m2, int m3):n2_(n2),n3_(n3),m2_(m2),m3_(m3),n1max_(n2+n3),n1min_(std::max(std::abs(n2-n3),std::abs(m2+m3))){
    
    ///possibly shift n1min_ if it is not aligned on even combinations of n1min+n2+n3
    n1min_+=(n2_+n3_+n1min_)%2;
    sz_=(int((n1max_-n1min_)/2)+1);
    data=std::vector<ftype>(sz_,0.0);
    ///Only order resulting in non-zero combinations
    m1_=-m2_-m3_;
    
     ///get Wigner3j (zero orders)
    Wigner3j<ftype> w3j0=Wigner3j<ftype>(n2_,n3_,0,0);    
    
    
    ///get Wigner3j coefficients 
    Wigner3j<ftype> w3jn=Wigner3j<ftype>(n2_,n3_,m2_,m3_);    
    size_t i=0; 
    for(int n = n1min_;n<=n1max_;n+=2){
        ftype scale=std::sqrt((2*n2_+1)*(2*n3_+1)/(4* M_PI)*(2*n+1));
        idx_[n]=i;///also fill out the index
        data[i++]=scale*w3j0[n]*w3jn[n];
    }



}

/*
 *  Reference for the Real Gaunt coefficients: Homeier, H.H.H., Steinborn, E.O., 1996. Some properties of the coupling coefficients of real spherical harmonics and their relation to Gaunt coefficients. Journal of Molecular Structure: THEOCHEM, Proceedings of the Second Electronic Computational Chemistry Conference 368, 31–37. https://doi.org/10.1016/S0166-1280(96)90531-X

 */
///@brief class to compute Real Gaunt coefficients (integral of 3 spherical harmonic base functions)
template<class ftype>
class GauntReal{
    public:
        GauntReal(){}
        GauntReal(int n2,int n3,int mu2,int mu3);
        std::vector<ftype> get ()const{return data;}
        ftype operator[](int n)const{return data[idx_.at(n)];}
        int nmin()const{return n1min_;}
        int nmax()const{return n1max_;}
        //int m()const{return m1_;}
    private:
    /* Unitary operator elements, see Homeier & Steinborn 1996 */
        static std::complex<ftype> Uel(int mu,int n,int m){
            if (mu ==  0){
                return kronecker(m,0)*kronecker(mu,0) + 0j;
            }else if (mu > 0){
                ///real elements only
                return M_SQRT1_2*(kronecker(m,mu)+csphase(m)*kronecker(m,-mu))+0j; 
            }else{
                ///mu < 0
                ///purely imaginary elements only
                return  1j*M_SQRT1_2*(csphase(m)*kronecker(m,mu)+kronecker(m,-mu)); 
            }
        }
        enum ABCcase{A,B,Bt,C}mcase;
        int n1max_=-1;
        int n1min_=-1;
        int mu1_[4];
        int n2_=-1;
        int mu2_=-1;
        int n3_=-1;
        int mu3_=-1;
        size_t sz_=0;
        std::vector<ftype> data;
        std::map<int,size_t> idx_;
};


template<class ftype>
GauntReal<ftype>::GauntReal(int n2,int n3,int mu2,int mu3):n2_(n2),n3_(n3),mu2_(mu2),mu3_(mu3){
    n1max_=n2+n2;

    n1min_=std::max(std::abs(n2-n3),std::min(std::abs(mu2+mu3),std::abs(mu2-mu3)));
    ///possibly shift by one to make sure even combinations of n1+n2+n3 are produced
    n1min_+=(n1max_+n1min_)%2;
    
    sz_=(n1max_-n1min_+1);
    ///obtain all non-zero cases for m1
    if(mu2 ==0 and mu3 ==0){
        mu1_[0]=0;
        mcase=ABCcase::C;
    }else if(mu3==0){
        mu1_[0]=-mu2;
        mu1_[1]=mu2;
        mcase=ABCcase::B;
    }else if(mu2==0){
        mu1_[0]=-mu3;
        mu1_[1]=mu3;

        mcase=ABCcase::Bt;
    }else{
        mu1_[0]=mu2+mu3;
        mu1_[1]=mu2-mu3;
        mu1_[2]=mu3-mu2;
        mu1_[3]=-mu3-mu2;
        mcase=ABCcase::A;
    }

    switch(mcase){
        case ABCcase::A:
                    //two valid contributions
                    Gaunt<ftype> gnt1(n2_,n3_,m2_,m3_);
                    Gaunt<ftype> gnt2(n2_,n3_,m2_,-m3_);

                     break;
        case ABCcase::B:
                    Gaunt<ftype> gnt(n2_,n3_,m2_,0);

                    break;
        case ABCcase::Bt:
                    std::vector<ftype)> gnt = Gaunt<ftype>(n2_,n3_,0,mu3_).get();

                    break;
        case ABCcase::C:
                    /// all orders are zero
                    data=Gaunt<ftype>(n2_,n3_,0,0).get();
                    break;
    
    }

}

            
#endif 
