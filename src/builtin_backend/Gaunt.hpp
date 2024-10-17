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
        ftype operator[](size_t i)const{return data[i];}
        int nmin()const{return n1min_;}
        int nmax()const{return n1max_;}
        std::vector<std::pair<int,int>> nm()const{return nm_;}
        std::vector<int> m()const{return mu1_;}
        inline size_t size() const{return data.size();}
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
                return M_SQRT1_2*(kronecker(m,mu)+csphase(m)*kronecker(m,-mu))+0i; 
                //according to sympy(https://docs.sympy.org/latest/modules/physics/wigner.html#sympy.physics.wigner.gaunt)
                //return M_SQRT1_2*(kronecker(m,mu)+csphase(mu)*kronecker(m,-mu))+0i; 
            }else{
                ///mu < 0
                ///purely imaginary elements only
                //according to Homeier et al 1996
                return  M_SQRT1_2*(csphase(m)*kronecker(m,mu)-kronecker(m,-mu))*1i; 
                //according to sympy(https://docs.sympy.org/latest/modules/physics/wigner.html#sympy.physics.wigner.gaunt)
                //return  M_SQRT1_2*(csphase(m)*kronecker(-m,mu)-kronecker(m,mu)*csphase(mu-m))*1i; 
            }
        }
        static inline int m_odd_factor(int m1,int m2,int m3){
        int modd=std::signbit(m1)+std::signbit(m2)+std::signbit(m3);
            if (modd%2 == 0){
                return 1;
            }else{
                return -1;
            }
        }
        int n1max_=-1;
        int n1min_=-1;
        std::vector<int> mu1_{};
        int n2_=-1;
        int n3_=-1;
        int mu2_=-1;
        int mu3_=-1;
        size_t sz_=0;
        std::vector<ftype> data;
        std::vector<std::pair<int,int>> nm_;/// Holds non-zero n1,m1 combinations
        std::map<std::pair<int,int>,size_t> idx_; ///For reverse lookups
};


template<class ftype>
GauntReal<ftype>::GauntReal(int n2,int n3,int mu2,int mu3):n2_(n2),n3_(n3),mu2_(mu2),mu3_(mu3),data(),nm_(){
    size_t idx=0;
   
    if (mu2_==0 && mu3_==0){
        //Eq 33 of Homeier
        //special case (quick return)
        Gaunt<ftype> gnt0(n2_,n3_,0,0);
        
        nm_.resize(gnt0.size(), std::make_pair(-1,-1));
        data.resize(gnt0.size(), 0.0);
        for (int n=gnt0.nmin();n<=gnt0.nmax();n+=2){
            //std::cout << n << " "<< m<< std::endl;
            nm_[idx]=std::make_pair(n,0);
            data[idx++]=gnt0[n];
        }
        //set minimum and maximum degree of valid Real Gaunt coefficients
        n1min_=gnt0.nmin();
        n1max_=gnt0.nmax();
        mu1_.push_back(0);   
        //no need to compute further
        return;

    }
    
    //create a list of non-zero options for mu1
    int mu1oa=mu2_+mu3_;
    if (m_odd_factor(mu1oa, mu2_, mu3_)==1){
        mu1_.push_back(mu1oa);
    }else{
        mu1_.push_back(-mu1oa);
    }

    if (mu2_*mu3_ != 0){
        //aditional options can be present
        int mu1ob=mu2_-mu3_;
        if (std::abs(mu1ob) != std::abs(mu1oa)){

            if (m_odd_factor(mu1ob, mu2_, mu3_)== 1){
                mu1_.push_back(mu1ob);
            }else{
                mu1_.push_back(-mu1ob);
            }
        }
    }
    

    Gaunt<ftype> gnt;
    ftype Un=0;

    ///reset min,max degree values (will bedetermined below)
    n1min_=std::numeric_limits<int>::max();
    n1max_=std::numeric_limits<int>::min();

    for (auto mu1 :mu1_){
        //std::cout << "C++ out "<< mu1 << " "<< mu2_ <<" " <<mu3_<< std::endl;

            //Try out the First term eq 31 Homeier (also includes the case for eq 32)
            /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
            Un= 2*csphase(mu1)*std::real(std::conj(Unitary(mu1,mu2_+mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,mu3_)));
            ///Un= 2*std::real(Unitary(mu1,mu2_+mu3_)*(Unitary(mu2_,mu2_)*Unitary(mu3_,mu3_)));
            if (Un != 0){
                gnt=Gaunt<ftype>(n2_,n3_,mu2_,mu3_);
            }else{
        
                //Try out the Second term eq 31 Homeier
                /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
                Un= 2*csphase(mu1)*std::real(std::conj(Unitary(mu1,mu2_-mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,-mu3_)));
                //Un= 2*std::real(Unitary(mu1,mu2_-mu3_)*(Unitary(mu2_,mu2_)*Unitary(mu3_,-mu3_)));
                gnt=Gaunt<ftype>(n2_,n3_,mu2_,-mu3_);

            }

    //std::cout << "C++ out "<< n1min_ << " "<< n1max_<< " " << mu1 << std::endl;

        //add the data and n,m pairs
        nm_.resize(gnt.size()+nm_.size(), std::make_pair(-1,-1));
        data.resize(gnt.size()+data.size(), 0.0);
        for (int n=gnt.nmin();n<=gnt.nmax();n+=2){
            //std::cout << n << " "<< m<< std::endl;
            nm_[idx]=std::make_pair(n,mu1);
            data[idx++]=Un*gnt[n];
        }
        //std::cout << "C++ out "<< gnt.nmin() << " "<< gnt.nmax()<< std::endl;
        //update min,max degrees present
        n1min_=std::min(n1min_,gnt.nmin());
        n1max_=std::max(n1max_,gnt.nmax());

    }
}
     

#endif 
