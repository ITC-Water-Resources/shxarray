/*! \file
 \brief Computation of (associated) Legendre functions (fast version)
 \copyright Roelof Rietbroek 2022
 \license
 This file is part of shxarray.
 */

#include <assert.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Legendre_nm.hpp"
#include <map>
#ifndef YNM_CPP_HPP
#define YNM_CPP_HPP

struct mni{
    int m;
    int n;
    size_t i;

};

//using mni=struct mni_str;

bool operator <(const mni& x, const mni& y) {
    return std::tie(x.m, x.n) < std::tie(y.m, y.n);
}

template <class ftype>
class Ynm_cpp{
   public:
    Ynm_cpp(const int nmax);
    Ynm_cpp(const size_t size, const int n [],const int m[ ]);
    Ynm_cpp(){};
    void set(const ftype lon, const ftype lat);
    

    //using mn=std::pair<int,int>;
    //inline ssize_t idx(int n, int m,int trig)const {
	//assert(m<=n);
	//int sgn=trig?-1:1;
	//mni mnisearch={m,n,
	//std::find(mnd.begin(), vec.end(), item) != vec.end()

	//return mnidx_.at(std::make_pair(sgn*m,n));
	
    //}

    int nmax() const { return legnm.nmax(); }
    int size() const { return sz_; }
    ftype * data(){return ynmdata_.data();}
    const ftype * data()const{return ynmdata_.data();}
    ftype & operator [](size_t i) {return ynmdata_[i];}
    const ftype & operator [](size_t i) const {return ynmdata_[i];}
    //inline std::map<mn,ssize_t> getmn()const{return mnidx_;}
    inline std::vector<mni> getmn()const{return mnidx_;}
   private:
    Legendre_nm<ftype> legnm;
    size_t sz_ = 0;
    std::vector<ftype> pnmcache_ = {};
    std::vector<ftype> ynmdata_ = {};
    std::vector<mni> mnidx_={};
    //std::map<mn,ssize_t> mnidx_={};
    ftype latprev=-1000; //initialize to impossible value
    bool sort=false;
};

template <class ftype>
Ynm_cpp<ftype>::Ynm_cpp(int nmax):legnm(nmax),sz_(2*(legnm.idx(nmax,nmax)+1)-(nmax+1))
		,pnmcache_(legnm.size()),ynmdata_(sz_,0.0),mnidx_(sz_){
/// Fill internal index	
    size_t i=0;
    for (int m=-nmax;m<=nmax;++m){
	for (int n=abs(m);n<=nmax;++n){
	    mnidx_[i]={m,n,i};
	    i++;
	}
    }
}

template<class ftype>
Ynm_cpp<ftype>::Ynm_cpp(const size_t size, const int n [],const int m[ ]):sz_(size),ynmdata_(size,0.0),mnidx_(sz_){
    
    ///find nmax
    int nmax=-1;
    //create a temporary map used to figure out the correct indices of the desired output order
    //map<std::pair<int,int>,ssize_t> mninputidx;
    for (size_t i=0;i<size;++i){
	nmax=(nmax<n[i])?n[i]:nmax;
	mnidx_[i]={m[i],n[i],i};
    }

    legnm=Legendre_nm<ftype>(nmax);
    pnmcache_=std::vector<ftype>(legnm.size());
    
    ///sort the vector so that order varies slowest (avoids unneeded recomputation of trigonometric functions)
    std::sort(mnidx_.begin(),mnidx_.end());

}



template <class ftype>
    void Ynm_cpp<ftype>::set(const ftype lon, const ftype lat){
    if (lat != latprev){
	///We need to recompute the associated legendre functions
	ftype costheta=sin(lat*M_PI/180.0);
	legnm.set(costheta,pnmcache_.data());
	//std::cout << "recompute for lat "<< lat << " " << latprev << std::endl;
	latprev=lat;
    }
    ftype lonr = lon*M_PI/180.0;
    ftype trig_mlon=0.0;
    int n,m;
    size_t idx;
    ///assign max to mold so it will trigger the computation of the trigonometric term on the first entry
    int mold=std::numeric_limits<int>::max();

    for (const auto & mni_item : mnidx_){
	m=mni_item.m;
	n=mni_item.n;
	idx=mni_item.i;
	///Note: the order (m) of the key varies slowest in the nmidx_ map so trigonometric factors are only recomputed when necessary 
	if (m != mold){
	    trig_mlon=(m<0)?sin(abs(m)*lonr):cos(m*lonr);
	    mold=m;
	}
	ynmdata_[idx]=trig_mlon*pnmcache_[legnm.idx(n,abs(m))];
    }


    }


#endif	/// YNM_CPP_HPP///
