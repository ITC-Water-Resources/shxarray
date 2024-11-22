/*! \file
 \brief 
 \copyright Roelof Rietbroek 2023
 \license
 This file is part of shxarray.
 */


#include <cmath>
#include <unordered_map>
#ifndef HELPERS_HPP_
#define HELPERS_HPP_


inline int csphase(int m){
    return 1-2*(std::abs(m)%2);
}

inline int kronecker(int n1,int n2){
   return n1==n2?1:0;
}




typedef std::pair<int, int> nmpair;
struct nm_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
         //we take advantage of the fact that we don't expect very large degrees which saturate the 32 bit integer
         //so we left shift the degree and add the absolute of the order (left shifted by 1 position) and add the signbit of the order
         // this should be fine for hasing degree order combinations up to 2**15 (32768)
        return static_cast<size_t>((pair.first<<15) + (std::abs(pair.second)<<1)+std::signbit(pair.second));
    }
};

typedef std::unordered_map<nmpair,size_t,nm_hash> nm_umap;

class Nmindex{
    public:
        Nmindex():nmmap_(1){};
        Nmindex(size_t maxsize):nmmap_(maxsize){};
        //size_t operator()(int n,int m)const{
            //return nmmap_[std::make_pair(n,m)];
        //}
        //const size_t & operator[](const nmpair & nm)const{
            //return nmmap_[nm];
        //}
        const size_t & operator[](const nmpair & nm)const{
            return nmmap_.at(nm);
        }
        void set(nmpair nm,size_t ix){
            nmmap_[nm]=ix;
        }
        //size_t get(const nmpair & nm){
            //return nmmap_[nm];
        //}
    private:
        nm_umap nmmap_;
};



#endif 
