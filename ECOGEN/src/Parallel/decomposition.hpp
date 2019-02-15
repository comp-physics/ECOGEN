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
    :nCells_(_nCells)
    {
        init( nProcs );
    }

    void init(int nProcs) noexcept
    {

        auto nCells_t=std::accumulate(nCells_.begin(),nCells_.end(),1,
                std::multiplies<int>());

        for(std::size_t d=0;d<Dim;++d)
        {
            const int level= static_cast<int >(std::log2(nCells_[d]))+1;
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
            processor_starts_per_level_.emplace_back(
                    std::vector<std::pair<key_type,int>>(1,
                        std::make_pair(key, nlocal)));

            while(count < nlocal)
            {
                if(is_valid(key)) ++count;
                ++key;
            }
        }
    }

    std::vector<key_type> get_keys( int _rank ) const noexcept
    {
        auto key=processor_starts_per_level_[_rank][0].first;
        auto nCells=processor_starts_per_level_[_rank][0].second;
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

    bool is_valid( const key_type& _key ) const noexcept
    {
        return (_key.coordinate()[0] < nCells_[0] &&
                _key.coordinate()[1] < nCells_[1] &&
                _key.coordinate()[2] < nCells_[2] ) ;
    }


    int get_load(const key_type& _k) const noexcept
    {
        return 1<< (_k.level()-base_level_);
    }


private:
    int base_level_=-1;
    std::array<int,Dim> nCells_;
    std::vector<std::vector<std::pair<key_type,int>>> processor_starts_per_level_;

};
}

#endif
