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
#ifndef FR_SH_WIGNER3j_HPP_
#define FR_SH_WIGNER3J_HPP_



/*
 * Reference for Wigner 3j recursion: Luscombe, J.H., Luban, M., 1998. Simplified recursive algorithm for Wigner 3j and 6j symbols. Physical Review E 57, 7274.
 * */
///@brief a class which computes Wigner3j symbols
template<class ftype>
    class Wigner3j{
        public:
            Wigner3j(){}
            Wigner3j(int j2,int j3,int m2, int m3);
	    const std::vector<ftype> get()const{return w3j_;}
            ftype operator[](int j)const{return w3j_[j-jmin_];}
        private:
            int j2_=-1;
            int j3_=-1;
            int m1_=-1;
            int m2_=-1;
            int m3_=-1;
            int jmin_=-1;
            int jmax_=-1;
            int sz_=0;
            std::vector<ftype> w3j_{};
            inline int jindex(int j)const{return j-jmin_;}
            void recursion();
            int signjmax()const{return csphase(j2_-j3_-m1_);}
            struct scales{
                ftype frac1;
                ftype frac2;
            };
           struct scales downward_scales(int j)const; 
           struct scales upward_scales(int j)const; 
    };



template<class ftype>
Wigner3j<ftype>::Wigner3j(int j2,int j3,int m2, int m3):j2_(j2),j3_(j3),m1_(-m2-m3),m2_(m2),m3_(m3),jmin_(std::max(std::abs(j2-j3),std::abs(m2+m3))),jmax_(j2+j3),sz_(jmax_-jmin_+1),w3j_(sz_,0.0)
{
    
    //std::cout << jmin_ <<" "<< jmax_ << " "<< sz_ << std::endl; 
    assert(jmax_ >= jmin_); 

    //Note vector is initialized to zero, so check here for quick returns of zeros
    //
    if (std::abs(m2_) > j2_ || std::abs(m3_) > j3_){
        return;
    }



    if (sz_ == 1){
    //This quick as only one term needs to be computed

        w3j_[0] = signjmax()/std::sqrt(2.0*jmin_+1.0); 
        return;
    }

    recursion();

}


template <class ftype>
void Wigner3j<ftype>::recursion(){
    ftype ratio=1.0;
    /* initialize values  at the end*/
    w3j_[jindex(jmax_)]=1e-4; /// !values at the boundaries are generally small !the scale 1d-3 is just an initial order of magnitude
    /* start with a downward recursion in the non-classical region*/
    int j; 
    for(j=jmax_;j!=jmin_+1;j--){
        
       struct scales factor= downward_scales(j);
        ratio=-factor.frac1/ratio-factor.frac2;
        //std::cout << "downward non-classical ratio " <<ratio<< std::endl; 
        if(std::abs(ratio) <= 1.0){
            /* local extrema reached, stop recursion*/
            
            break;
        }
        w3j_[jindex(j-1)]=ratio*w3j_[jindex(j)];
    }
    
    int jupper_classical=j;
    ftype ValueUpper_classical=w3j_[jindex(jupper_classical)];

    /* upward recursion in non-classical region */


    ///initial bottom value (copy value from the top)
    ratio=1.0;
    w3j_[jindex(jmin_)]=w3j_[jindex(jmax_)];
    for(j=jmin_;j!=jupper_classical;j++){
        
       struct scales factor= upward_scales(j);
        ratio=-factor.frac1/ratio-factor.frac2;
        //std::cout << "upward non-classical ratio" <<ratio<< std::endl; 
        if(std::abs(ratio) < 1.0){
            /* local extrema reached, stop recursion*/
            
            break;
        }
        w3j_[jindex(j+1)]=ratio*w3j_[jindex(j)];
    }
    
    int jlower_classical=j;

    ftype ValueLower_classical=w3j_[jindex(jlower_classical)];

    /* perform the classical recursion in the direction which has the least amount of total recursion steps*/
    if (jmax_-jlower_classical < jupper_classical-jmin_){
    /*use downward classical recursion */   
        for(j=jupper_classical;j!=jlower_classical;j--){
            
            struct scales factor= downward_scales(j);
            //std::cout << "downward classical recursion " <<j<<" "<<factor.frac1 <<" "<<factor.frac2<< std::endl; 
            w3j_[jindex(j-1)]=-factor.frac2*w3j_[jindex(j)]-factor.frac1*w3j_[jindex(j+1)];
        }

        /* rescale lower part to match scale at jlower_classical*/
        ftype scale=w3j_[jindex(jlower_classical)]/ValueLower_classical;
        for(int i=0;i<jindex(jlower_classical);i++){
            w3j_[i]*=scale;
        }

    }else{
    /*use upward classical recursion */   

        for(j=jlower_classical;j<jupper_classical;j++){
            
        //std::cout << "upward classical recursion" <<j<< std::endl; 
           struct scales factor= upward_scales(j);
            w3j_[jindex(j+1)]=-factor.frac2*w3j_[jindex(j)]-factor.frac1*w3j_[jindex(j-1)];
        }
    

        /* rescale upper part to match scale at jupper_classical*/
        ftype scale=w3j_[jindex(jupper_classical)]/ValueUpper_classical;
        for(int i=jindex(jupper_classical+1);i<=jindex(jmax_);i++){
            w3j_[i]*=scale;
        }
    }
    /* apply an additional scaling to normalize the whole vector */

    ftype norm = 0.0;
    for (int j=jmin_;j<=jmax_;j++){
        norm+=static_cast<ftype>(2*j+1)*std::pow(w3j_[jindex(j)],2);
    }
    norm=1.0/std::sqrt(norm);
    ///Also make sure the sign at the end fits 
    //std::cout << signjmax() << " norm " << w3j_[jindex(jmax_)] << std::endl;
    if ( (signjmax()*w3j_[jindex(jmax_)] ) < 0){
        //adapt the sign so that the last value will have the correct sign
        norm*=-1.0;
        //std::cout << "negating norm " << norm << std::endl;
    }

    for(auto && it: w3j_){
        it*=norm;
    }

}


template <class ftype>
struct Wigner3j<ftype>::scales Wigner3j<ftype>::downward_scales(int j)const{
    struct scales x_y_over_z;
    /*common denominator part*/
    ftype tmp1=j*j-std::pow(j2_-j3_,2);
    tmp1*=(std::pow(j2_+j3_+1,2)-j*j);
    tmp1*=(j*j-m1_*m1_);
    //return x_y_over_z;
    /*xratio*/
    ftype tmp2=(std::pow(j+1,2)-std::pow(j2_-j3_,2));
    tmp2*=(std::pow(j2_+j3_+1,2)-std::pow(j+1,2));
    tmp2*=(std::pow(j+1,2)-(m1_*m1_));
    
    /* compute x_over_z and put in frac1*/
    x_y_over_z.frac1=(static_cast<ftype>(j)/(j+1))*std::sqrt(tmp2/tmp1);///x_over_z

    /*y over z ratio*/

   
    tmp2=((-m1_)*(j2_*(j2_+1)-j3_*(j3_+1))-(m2_-m3_)*(j*(j+1)))/std::sqrt(tmp1);

    x_y_over_z.frac2=(1.0+static_cast<ftype>(j)/(j+1))*tmp2;///y_over_z

    return x_y_over_z;

}


template <class ftype>
struct Wigner3j<ftype>::scales Wigner3j<ftype>::upward_scales(int j)const{
    struct scales z_y_over_x;
    /* quick return */
    if(j == 0){ 
       z_y_over_x.frac1=0.0; ///z_over_x
       z_y_over_x.frac2=-(m2_-m3_)/std::sqrt(std::pow(j2_+j3_+1,2)-1.0);///y_over_x
       return z_y_over_x;
    }

    /*common denominator part*/
    ftype tmp1=std::pow(j+1,2) - std::pow(j2_-j3_,2);
    tmp1*=(std::pow(j2_+j3_+1,2) -std::pow(j+1,2));
    tmp1*=(std::pow(j+1,2) - m1_*m1_);

    /* zratio */
    ftype tmp2=(j*j-std::pow(j2_-j3_,2));
    tmp2*=(std::pow(j2_+j3_+1,2)-(j*j));
    tmp2*=((j*j)-(m1_*m1_));

    z_y_over_x.frac1=(static_cast<ftype>(j+1)/(j))*std::sqrt(tmp2/tmp1);
    
    /*yratio */
    tmp2=((-m1_)*(j2_*(j2_+1)-j3_*(j3_+1))-(m2_-m3_)*(j*(j+1)))/std::sqrt(tmp1);

    z_y_over_x.frac2=(2.0+1.0/j)*tmp2;

    return z_y_over_x;

}



            
#endif 
