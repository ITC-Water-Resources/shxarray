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
#include <unordered_map>
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
        Gaunt(int nmax);
        Gaunt(int n2,int n3,int m2,int m3);
        std::vector<ftype> get ()const{return data;}
        void set(const int n2,const int n3,const int m2,const int m3);
        ftype operator[](int n)const{return data[idx(n)];}
        int nmin()const{return n1min_;}
        int nmax()const{return n1max_;}
        int m()const{return m1_;}
        inline size_t size() const{return sz_;}
        inline size_t idx(int n)const{
            assert(n >0);
            assert (n >= n1min_);
            assert (n <= n1max_);
            assert ((n-n1min_)%2 ==0);
            return (n-n1min_)/2;}
    private:
        int m1_=-1;
        int n2_=-1;
        int n3_=-1;
        int m2_=-1;
        int m3_=-1;
        int n1max_=-1;
        int n1min_=-1;
        size_t sz_=0;
        Wigner3j<ftype> w3j0_;
        Wigner3j<ftype> w3jn_;
        std::vector<ftype> data;
        //std::unordered_map<int,size_t> idx_;
};



template<class ftype>
Gaunt<ftype>::Gaunt(int nmax):n1max_(nmax),n1min_(0),w3j0_(nmax),w3jn_(nmax),data(n1max_+1){

    
}

template<class ftype>
Gaunt<ftype>::Gaunt(int n2,int n3, int m2, int m3):n2_(n2),n3_(n3),m2_(m2),m3_(m3),n1max_(n2+n3),n1min_(std::max(std::abs(n2-n3),std::abs(m2+m3))),w3j0_(n1max_),w3jn_(n1max_),data(n1max_+1){
    set(n2,n3,m2,m3);
}

template<class ftype>
void Gaunt<ftype>::set(const int n2,const int n3, const int m2, const int m3){
    n2_=n2;
    n3_=n3;
    m2_=m2;
    m3_=m3;
    n1max_=(n2+n3);
    n1min_=(std::max(std::abs(n2-n3),std::abs(m2+m3)));
    
    ///possibly shift n1min_ if it is not aligned on even combinations of n1min+n2+n3
    n1min_+=(n2_+n3_+n1min_)%2;
    sz_=(int((n1max_-n1min_)/2)+1);
    data.resize(sz_,0.0);
    ///Only order resulting in non-zero combinations
    //m1_=m2_+m3_;
    m1_=-m2_-m3_;
    
     ///get Wigner3j (zero orders)
    w3j0_.set(n2_,n3_,0,0);    
    
    
    ///get Wigner3j coefficients 
    w3jn_.set(n2_,n3_,m2_,m3_);    
    
    for(int n = n1min_;n<=n1max_;n+=2){
        ftype scale=std::sqrt((2*n2_+1)*(2*n3_+1)/(4* M_PI)*(2*n+1));
        data[idx(n)]=scale*w3j0_[n]*w3jn_[n];
    }



}


/*
*  Reference for the Real Gaunt coefficients: Homeier, H.H.H., Steinborn, E.O., 1996. Some properties of the coupling coefficients of real spherical harmonics and their relation to Gaunt coefficients. Journal of Molecular Structure: THEOCHEM, Proceedings of the Second Electronic Computational Chemistry Conference 368, 31–37. https://doi.org/10.1016/S0166-1280(96)90531-X

 */
///@brief class to compute Real Gaunt coefficients (integral of 3 real spherical harmonic base functions)
template<class ftype>
class GauntReal{
    public:
        GauntReal();
        GauntReal(int nmax);
        GauntReal(int n2,int n3,int mu2,int mu3);
        std::vector<ftype> get ()const{return data;}
        void set(const int n2,const int n3,const int mu2,const int mu3);
        ftype operator[](size_t i)const{return data[i];}
        ftype at(int n1,int m1)const{
            return data[idx(n1,m1)];}
        int nmin()const{return std::min(n1min_[0],n1min_[1]);}
        int nmax()const{return std::max(n1max_[0],n1max_[1]);}
        const std::vector<nmpair> & nmvec()const{return nm_;}
        const nmpair & nm(size_t i)const {return nm_[i];}
        //std::vector<int> m()const{return mu1_;}
        inline size_t size() const{return data.size();}
        inline size_t idx(int n,int m)const{
            
            size_t im=100; //invalid on purpose
            size_t shft=0;
            if (m == mu1_[0]){
                im=0;
            }else if (m == mu1_[1]){
                im=1;
                shft=gnt_[0].size();
            }
            assert(im <2);
            assert(n> 0 and n>= n1min_[im] and n <= n1max_[im]);

            return shft+gnt_[im].idx(n); 
        }

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
        int n1max_[2]={std::numeric_limits<int>::min(),std::numeric_limits<int>::min()};
        int n1min_[2]={std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
        int mu1_[2]={-1,-1};
        size_t ithm=0;
        int n2_=-1;
        int n3_=-1;
        int mu2_=-1;
        int mu3_=-1;
        Gaunt<ftype> gnt_[2];
        std::vector<ftype> data;
        std::vector<nmpair> nm_;/// Holds non-zero n1,m1 combinations
};

template<class ftype>
GauntReal<ftype>::GauntReal(){
}

template<class ftype>
GauntReal<ftype>::GauntReal(int nmax):data((nmax+1)*2),nm_((nmax+1)*2){
    /// Constructor to reserve(allocate) for maximum capacity
}
template<class ftype>
GauntReal<ftype>::GauntReal(int n2,int n3,int mu2,int mu3):n2_(n2),n3_(n3),mu2_(mu2),mu3_(mu3),data(),nm_(){
    set(n2,n3,mu2,mu3);
}

template<class ftype>
void GauntReal<ftype>::set(const int n2,const int n3,const int mu2,const int mu3){
    n2_=n2;
    n3_=n3;
    mu2_=mu2;
    mu3_=mu3;
    data.resize(0);
    nm_.resize(0);
    ithm=0;
    
    ///reset min,max degree values (will bedetermined below)
    /// this reset is needed for functioning nmin() and nmax() functions
    for (size_t i=0;i<2;++i){
        n1max_[i]=std::numeric_limits<int>::min();
        n1min_[i]=std::numeric_limits<int>::max();
    }
   
    if (mu2_==0 && mu3_==0){
    
        //Eq 33 of Homeier
        //special case (quick return)
        gnt_[0].set(n2_,n3_,0,0);
        
        nm_.resize(gnt_[0].size(), std::make_pair(-1,-1));
        data.resize(gnt_[0].size(), 0.0);
        //set minimum and maximum degree of valid Real Gaunt coefficients
        mu1_[ithm]=0;
        n1min_[ithm]=gnt_[0].nmin();
        n1max_[ithm]=gnt_[0].nmax();
        for (int n=n1min_[ithm];n<=n1max_[ithm];n+=2){
            //std::cout << n << " "<< m<< std::endl;
            size_t ix=idx(n,0);
            nm_[ix]=std::make_pair(n,0);
            data[ix]=gnt_[0][n];
            
        }

        //no need to compute further
        return;

    }
    
    //create a list of non-zero options for mu1
    int mu1oa=mu2_+mu3_;
    if (m_odd_factor(mu1oa, mu2_, mu3_)==1){
        mu1_[ithm]=mu1oa;
    }else{
        mu1_[ithm]=-mu1oa;
    }

    if (mu2_*mu3_ != 0){
        //aditional options can be present
        int mu1ob=mu2_-mu3_;
        if (std::abs(mu1ob) != std::abs(mu1oa)){
            ithm=1;
            if (m_odd_factor(mu1ob, mu2_, mu3_)== 1){
                mu1_[ithm]=mu1ob;
            }else{
                mu1_[ithm]=-mu1ob;
            }
        }
    }
    

    ftype Un=0;

    for (size_t im=0;im<=ithm;++im){
        //std::cout << "C++ out "<< mu1 << " "<< mu2_ <<" " <<mu3_<< std::endl;
            int mu1=mu1_[im];
            //Try out the First term eq 31 Homeier (also includes the case for eq 32)
            /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
            Un= 2*csphase(mu1)*std::real(std::conj(Unitary(mu1,mu2_+mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,mu3_)));
            ///Un= 2*std::real(Unitary(mu1,mu2_+mu3_)*(Unitary(mu2_,mu2_)*Unitary(mu3_,mu3_)));
            if (Un != 0){
                gnt_[im].set(n2_,n3_,mu2_,mu3_);
            }else{
        
                //Try out the Second term eq 31 Homeier
                /// Note because Homeier et al define Gaunt coefficients a bit differently we need to additionally apply the csphase 
                Un= 2*csphase(mu1)*std::real(std::conj(Unitary(mu1,mu2_-mu3_))*(Unitary(mu2_,mu2_)*Unitary(mu3_,-mu3_)));
                //Un= 2*std::real(Unitary(mu1,mu2_-mu3_)*(Unitary(mu2_,mu2_)*Unitary(mu3_,-mu3_)));
                gnt_[im].set(n2_,n3_,mu2_,-mu3_);

            }

    //std::cout << "C++ out "<< n1min_ << " "<< n1max_<< " " << mu1 << std::endl;

        //add the data and n,m pairs
        nm_.resize(gnt_[im].size()+nm_.size(), std::make_pair(-1,-1));
        data.resize(gnt_[im].size()+data.size(), 0.0);
        //update min,max degrees present
        n1min_[im]=gnt_[im].nmin();
        n1max_[im]=gnt_[im].nmax();

        for (int n=n1min_[im];n<=n1max_[im];n+=2){
            //std::cout << n << " "<< m<< std::endl;
            size_t ix=idx(n,mu1);
            nm_[ix]=std::make_pair(n,mu1);
            data[ix]=Un*gnt_[im][n];
        }

    }
}
     

#endif 
