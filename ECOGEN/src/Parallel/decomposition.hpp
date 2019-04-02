//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDED_DECOMP_HPP
#define INCLUDED_DECOMP_HPP

//! \file      decomposition.hpp
//! \author    B. Dorschner
//! \version   1.0
//! \date      February 15 2019

#include <algorithm>
#include <numeric>
#include <functional>
#include "key.hpp"

namespace decomposition
{

class Decomposition
{

public:
    static constexpr int Dim=3;
    using key_type=Key<Dim>;


public: //Ctors
	Decomposition() = default;
	Decomposition(const Decomposition& other) = default;
	Decomposition(Decomposition&& other) = default;
	Decomposition& operator=(const Decomposition& other) & = default;
	Decomposition& operator=(Decomposition&& other) & = default;
	~Decomposition() = default;



    Decomposition(std::array<int,Dim> _nCells, int nProcs)
    :nCells_global_(_nCells)
    {
        init( nProcs );
    }

    void init(int nProcs) noexcept
    {

        auto nCells_t=std::accumulate(nCells_global_.begin(),nCells_global_.end(),1,
                std::multiplies<int>());

        for(std::size_t d=0;d<Dim;++d)
        {
            const int level= static_cast<int >(std::log2(nCells_global_[d]))+1;
            if(level>base_level_) base_level_=level;
        }


        float chunks = static_cast<float>(nCells_t)/nProcs;
        key_type key(0,0,0,base_level_);
        for ( int i=0; i<nProcs;++i )
        {
            size_t start= (i*chunks);
            size_t end= ((i+1)*chunks);
            const int nlocal=end-start;
            int count=0;
            key_rank_map_.emplace_back(key);
            nCells_per_rank_.emplace_back(nlocal);

            while(count < nlocal)
            {
                if(is_valid(key)) ++count;
                ++key;
            }

            //Store also end, to check validity:
            if(i==nProcs-1)
            {
                nCells_per_rank_.emplace_back(0);
                key_rank_map_.emplace_back(key);
            }
        }
    }

    std::vector<key_type> get_keys( int _rank ) const noexcept
    {
        auto key=key_rank_map_[_rank];
        auto nCells=nCells_per_rank_[_rank];
         std::vector<key_type> keys(nCells);

        for(int i =0; i<nCells;++i)
        {
            bool valid=false;
            while(!valid)
            {
                if(is_valid(key))  
                {
                    valid=true;
                    keys[i]=key;
                }
                ++key;
            }
        }
        return keys;
    }

    int get_rank(const key_type _key)
    {
      for(std::size_t i=0; i< key_rank_map_.size();++i) 
      {
          if(_key < key_rank_map_[i] )
              return i-1;
      }
     return -1;
    }

    bool is_valid( const key_type& _key ) const noexcept
    {
        return (_key.coordinate()[0] < nCells_global_[0] &&
                _key.coordinate()[1] < nCells_global_[1] &&
                _key.coordinate()[2] < nCells_global_[2] ) ;
    }


    int get_load(const key_type& _k) const noexcept
    {
        return 1<< (_k.level()-base_level_);
    }

    
    int& base_level()noexcept{return base_level_;}
    const int& base_level()const noexcept{return base_level_;}

    template<class Coord>
    bool is_inside(const Coord& _coord)
    {
       for(int d=0;d<Dim;++d)
       {
           if( _coord[d]<0 || _coord[d]>= nCells_global_[d])
               return false;
       }
       return true;
    }


    const key_type& start_key(int i) const noexcept{return key_rank_map_[i];}
    key_type& start_key(int i) noexcept{return key_rank_map_[i];}


private:
    int base_level_=-1;
    std::array<int,Dim> nCells_global_;
    std::vector<key_type> key_rank_map_;
    std::vector<int> nCells_per_rank_; //Not updated once the starts are changed

};
}

#endif
