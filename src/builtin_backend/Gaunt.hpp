/*! \file
 \brief 
 \copyright Roelof Rietbroek 2024
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
#include "Wigner3j.hpp"

#ifndef FR_SH_GAUNT_HPP_
#define FR_SH_GAUNT_HPP_


///@brief class to compute Gaunt coefficients (integral of 3 (complex) spherical harmonic base functions)
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
        inline size_t size() const{return sz_;}
    private:
        int m1_=-1;
        int n2_=-1;
        int n3_=-1;
        int m2_=-1;
        int m3_=-1;
        int n1max_=-1;
        int n1min_=-1;
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
    //m1_=m2_+m3_;
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
///@brief class to compute Real Gaunt coefficients (integral of 3 real spherical harmonic base functions)
template<class ftype>
class GauntReal{
    public:
        GauntReal(){}
        GauntReal(int n2,int n3,int mu2,int mu3);
        std::vector<ftype> get ()const{return data;}
        //ftype operator[](int n)const{return data[idx_.at(n)];}
        int nmin()const{return n1min_;}
        int nmax()const{return n1max_;}
        std::vector<std::pair<int,int>> nm()const{return nm_;}
        int m()const{return mu1_;}
    private:
        //[> Unitary operator elements, see Homeier & Steinborn 1996 <]
        static std::complex<ftype> Unitary(int mu,int m){
            using namespace std::complex_literals;
            if (std::abs(mu) != std::abs(m)){
                return 0.0+0i;
            }
            else if (mu ==  0){
                return 1.0 + 0i;
            }else if (mu > 0){
                ///real elements only
                //according to Homeier et al 1996
                //return M_SQRT1_2*(kronecker(m,mu)+csphase(m)*kronecker(m,-mu))+0j; 
                //according to sympy(https://docs.sympy.org/latest/modules/physics/wigner.html#sympy.physics.wigner.gaunt)
                return M_SQRT1_2*(kronecker(m,mu)+csphase(mu)*kronecker(m,-mu))+0i; 
            }else{
                ///mu < 0
                ///purely imaginary elements only
                //according to Homeier et al 1996
                //return  1j*M_SQRT1_2*(csphase(m)*kronecker(m,mu)+kronecker(m,-mu)); 
                //according to sympy(https://docs.sympy.org/latest/modules/physics/wigner.html#sympy.physics.wigner.gaunt)
                return  M_SQRT1_2*(csphase(m)*kronecker(-m,mu)-kronecker(m,mu)*csphase(mu-m))*1i; 
            }
        }
        enum ABCcase{A,B,Bt,C}mcase;
        int n1max_=-1;
        int n1min_=-1;
        int mu1_=-1;
        int n2_=-1;
        int n3_=-1;
        int mu2_=-1;
        int mu3_=-1;
        size_t sz_=0;
        std::vector<ftype> data;
        std::vector<std::pair<int,int>> nm_;/// Holds non-zero n1,m1 combinations
        //std::map<int,size_t> idx_;
};


template<class ftype>
GauntReal<ftype>::GauntReal(int n2,int n3,int mu2,int mu3):n2_(n2),n3_(n3),mu2_(mu2),mu3_(mu3),data(),nm_(){
    size_t idx=0;
   
    if (mu2_==0 && mu3_==0){
        //special case (quick return)
        Gaunt<ftype> gnt0(n2_,n3_,0,0);
        
        nm_.resize(gnt0.size(), std::make_pair(-1,-1));
        data.resize(gnt0.size(), 0.0);
        n1min_=gnt0.nmin();
        n1max_=gnt0.nmax();
        mu1_=0;
        for (int n=n1min_;n<=n1max_;n+=2){
            //std::cout << n << " "<< m<< std::endl;
            nm_[idx]=std::make_pair(n,0);
            data[idx++]=gnt0[n];
        }
        //no need to compute further
        return;

    }

    ///Try the first option for mu1
    int mu1v1 = -mu2_-mu3_;
    /// Zero values occur when an odd number of negative orders are present (so skip that case)
    if (!std::signbit(mu1v1*mu2_*mu3_)){
        
        Gaunt<ftype> gnt1(n2_,n3_,mu2_,mu3_);
        //std::cout << "computing mu1 v1 " << " " << mu1v1 << std::endl;
        n1min_=gnt1.nmin();
        n1max_=gnt1.nmax();
        mu1_=mu1v1;
        /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
        ftype Un1= 2*csphase(mu1v1)*std::real(std::conj(Unitary(mu1v1,mu2_+mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,mu3_)));
        
        
        nm_.resize(gnt1.size(), std::make_pair(-1,-1));
        data.resize(gnt1.size(), 0.0);
        for (int n=gnt1.nmin();n<=gnt1.nmax();n+=2){
            //std::cout << n << " "<< m<< std::endl;
            nm_[idx]=std::make_pair(n,mu1v1);
            data[idx++]=Un1*gnt1[n];
        }
        //this mu1 option is mutually exclusive with the next one so no need to try the next option
        return;     
        
    }
    
    ///Second possibility for mu1
    int mu1v2= -mu2_ +mu3_;
    if (!std::signbit(mu1v2*mu2_*mu3_)){
        
        //std::cout << "computing mu1 v2 " << " " << mu1v2 << std::endl;
        Gaunt<ftype> gnt2(n2_,n3_,mu2_,-mu3_);
        
        n1min_=gnt2.nmin();
        n1max_=gnt2.nmax();
        mu1_=mu1v2;
        
        /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
        
        ftype Un2= 2*csphase(mu1v2)*std::real(std::conj(Unitary(mu1v2,mu2_-mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,-mu3_)));
        
        nm_.resize(gnt2.size()+nm_.size(), std::make_pair(-1,-1));
        data.resize(gnt2.size()+data.size(), 0.0);
        for (int n=gnt2.nmin();n<=gnt2.nmax();n+=2){
            //std::cout << n << " "<< m<< std::endl;
            nm_[idx]=std::make_pair(n,mu1v2);
            data[idx++]=Un2*gnt2[n];
        }

        
    }
    
}

#endif 
