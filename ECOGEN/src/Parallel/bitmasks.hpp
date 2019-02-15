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

#ifndef INCLUDED_BIT_MASKS_HPP
#define INCLUDED_BIT_MASKS_HPP

//! \file      bitmasks.hpp
//! \author    B. Dorschner
//! \version   1.0
//! \date      February 15 2019

#include <stdint.h>
#include <array>

namespace decomposition
{

namespace bitmask_t
{

using index_t = unsigned long long int;
using difference_type = long long int;
using scalar_coordinate_type = int;
using level_type = int;
constexpr level_type max_level = 19;

//Note 0b indicates that the binary literal is used
//The key bits are encoded as follows
//0--19*3-1=0--56 bit: coordinates
//57--61 bit: level
//61--63 bit: flags
//                                    --1--2--3--4--5--6--7--8--9-10-11-12-13-14-15-16-17-18-19--lev--
//                                    --|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|----|-|
 constexpr index_t level_mask    = 0b0000000000000000000000000000000000000000000000000000000001111100;
 constexpr index_t coord_mask    = 0b1111111111111111111111111111111111111111111111111111111110000000;
 constexpr index_t flag_mask     = 0b0000000000000000000000000000000000000000000000000000000000000011;
 constexpr index_t lo_mask       = 0b0000000000000000000000000000000000000000000000000000000000000001;
 constexpr index_t hi_mask       = 0b1000000000000000000000000000000000000000000000000000000000000000;
 constexpr index_t hi_3_mask     = 0b1110000000000000000000000000000000000000000000000000000000000000;
 constexpr index_t x_mask        = 0b0010010010010010010010010010010010010010010010010010010010000000;
 constexpr index_t y_mask        = 0b0100100100100100100100100100100100100100100100100100100100000000;
 constexpr index_t z_mask        = 0b1001001001001001001001001001001001001001001001001001001000000000;

 constexpr index_t min_0         = 0b0000000000000000000000000000000000000000000000000000000000000001;
 constexpr index_t min_1         = 0b0000000000000000000000000000000000000000000000000000000000000101;
 constexpr index_t min_2         = 0b0000000000000000000000000000000000000000000000000000000000001001;
 constexpr index_t min_3         = 0b0000000000000000000000000000000000000000000000000000000000001101;
 constexpr index_t min_4         = 0b0000000000000000000000000000000000000000000000000000000000010001;
 constexpr index_t min_5         = 0b0000000000000000000000000000000000000000000000000000000000010101;
 constexpr index_t min_6         = 0b0000000000000000000000000000000000000000000000000000000000011001;
 constexpr index_t min_7         = 0b0000000000000000000000000000000000000000000000000000000000011101;
 constexpr index_t min_8         = 0b0000000000000000000000000000000000000000000000000000000000100001;
 constexpr index_t min_9         = 0b0000000000000000000000000000000000000000000000000000000000100101;
 constexpr index_t min_10        = 0b0000000000000000000000000000000000000000000000000000000000101001;
 constexpr index_t min_11        = 0b0000000000000000000000000000000000000000000000000000000000101101;
 constexpr index_t min_12        = 0b0000000000000000000000000000000000000000000000000000000000110001;
 constexpr index_t min_13        = 0b0000000000000000000000000000000000000000000000000000000000110101;
 constexpr index_t min_14        = 0b0000000000000000000000000000000000000000000000000000000000111001;
 constexpr index_t min_15        = 0b0000000000000000000000000000000000000000000000000000000000111101;
 constexpr index_t min_16        = 0b0000000000000000000000000000000000000000000000000000000001000001;
 constexpr index_t min_17        = 0b0000000000000000000000000000000000000000000000000000000001000101;
 constexpr index_t min_18        = 0b0000000000000000000000000000000000000000000000000000000001001001;
 constexpr index_t min_19        = 0b0000000000000000000000000000000000000000000000000000000001001101;

 constexpr std::array<index_t,20> min_arr = {{ min_0,  min_1,  min_2,  min_3,  min_4,
	                                                 min_5,  min_6,  min_7,  min_8,  min_9,
	                                                 min_10, min_11, min_12, min_13, min_14,
	                                                 min_15, min_16, min_17, min_18, min_19  }};

 constexpr index_t max_0         = 0b0000000000000000000000000000000000000000000000000000000000000001;
 constexpr index_t max_1         = 0b1110000000000000000000000000000000000000000000000000000000000101;
 constexpr index_t max_2         = 0b1111110000000000000000000000000000000000000000000000000000001001;
 constexpr index_t max_3         = 0b1111111110000000000000000000000000000000000000000000000000001101;
 constexpr index_t max_4         = 0b1111111111110000000000000000000000000000000000000000000000010001;
 constexpr index_t max_5         = 0b1111111111111110000000000000000000000000000000000000000000010101;
 constexpr index_t max_6         = 0b1111111111111111110000000000000000000000000000000000000000011001;
 constexpr index_t max_7         = 0b1111111111111111111110000000000000000000000000000000000000011101;
 constexpr index_t max_8         = 0b1111111111111111111111110000000000000000000000000000000000100001;
 constexpr index_t max_9         = 0b1111111111111111111111111110000000000000000000000000000000100101;
 constexpr index_t max_10        = 0b1111111111111111111111111111110000000000000000000000000000101001;
 constexpr index_t max_11        = 0b1111111111111111111111111111111110000000000000000000000000101101;
 constexpr index_t max_12        = 0b1111111111111111111111111111111111110000000000000000000000110001;
 constexpr index_t max_13        = 0b1111111111111111111111111111111111111110000000000000000000110101;
 constexpr index_t max_14        = 0b1111111111111111111111111111111111111111110000000000000000111001;
 constexpr index_t max_15        = 0b1111111111111111111111111111111111111111111110000000000000111101;
 constexpr index_t max_16        = 0b1111111111111111111111111111111111111111111111110000000001000001;
 constexpr index_t max_17        = 0b1111111111111111111111111111111111111111111111111110000001000101;
 constexpr index_t max_18        = 0b1111111111111111111111111111111111111111111111111111110001001001;
 constexpr index_t max_19        = 0b1111111111111111111111111111111111111111111111111111111111001101;

 constexpr std::array<index_t,20> max_arr = {{ max_0,  max_1,  max_2,  max_3,  max_4,
	                                                 max_5,  max_6,  max_7,  max_8,  max_9,
	                                                 max_10, max_11, max_12, max_13, max_14,
	                                                 max_15, max_16, max_17, max_18, max_19  }};

 constexpr index_t coord_mask_0  = 0b0000000000000000000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_1  = 0b1110000000000000000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_2  = 0b1111110000000000000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_3  = 0b1111111110000000000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_4  = 0b1111111111110000000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_5  = 0b1111111111111110000000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_6  = 0b1111111111111111110000000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_7  = 0b1111111111111111111110000000000000000000000000000000000000000000;
 constexpr index_t coord_mask_8  = 0b1111111111111111111111110000000000000000000000000000000000000000;
 constexpr index_t coord_mask_9  = 0b1111111111111111111111111110000000000000000000000000000000000000;
 constexpr index_t coord_mask_10 = 0b1111111111111111111111111111110000000000000000000000000000000000;
 constexpr index_t coord_mask_11 = 0b1111111111111111111111111111111110000000000000000000000000000000;
 constexpr index_t coord_mask_12 = 0b1111111111111111111111111111111111110000000000000000000000000000;
 constexpr index_t coord_mask_13 = 0b1111111111111111111111111111111111111110000000000000000000000000;
 constexpr index_t coord_mask_14 = 0b1111111111111111111111111111111111111111110000000000000000000000;
 constexpr index_t coord_mask_15 = 0b1111111111111111111111111111111111111111111110000000000000000000;
 constexpr index_t coord_mask_16 = 0b1111111111111111111111111111111111111111111111110000000000000000;
 constexpr index_t coord_mask_17 = 0b1111111111111111111111111111111111111111111111111110000000000000;
 constexpr index_t coord_mask_18 = 0b1111111111111111111111111111111111111111111111111111110000000000;
 constexpr index_t coord_mask_19 = 0b1111111111111111111111111111111111111111111111111111111110000000;

 constexpr std::array<index_t,20> coord_mask_arr = {{ coord_mask_0,  coord_mask_1,  coord_mask_2,
                                                            coord_mask_3,  coord_mask_4,  coord_mask_5,
                                                            coord_mask_6,  coord_mask_7,  coord_mask_8,
                                                            coord_mask_9,  coord_mask_10, coord_mask_11,
                                                            coord_mask_12, coord_mask_13, coord_mask_14,
                                                            coord_mask_15, coord_mask_16, coord_mask_17,
                                                            coord_mask_18, coord_mask_19  }};

 constexpr scalar_coordinate_type max_coord_0  =      1;
 constexpr scalar_coordinate_type max_coord_1  =      2;
 constexpr scalar_coordinate_type max_coord_2  =      4;
 constexpr scalar_coordinate_type max_coord_3  =      8;
 constexpr scalar_coordinate_type max_coord_4  =     16;
 constexpr scalar_coordinate_type max_coord_5  =     32;
 constexpr scalar_coordinate_type max_coord_6  =     64;
 constexpr scalar_coordinate_type max_coord_7  =    128;
 constexpr scalar_coordinate_type max_coord_8  =    256;
 constexpr scalar_coordinate_type max_coord_9  =    512;
 constexpr scalar_coordinate_type max_coord_10 =   1024;
 constexpr scalar_coordinate_type max_coord_11 =   2048;
 constexpr scalar_coordinate_type max_coord_12 =   4096;
 constexpr scalar_coordinate_type max_coord_13 =   8192;
 constexpr scalar_coordinate_type max_coord_14 =  16384;
 constexpr scalar_coordinate_type max_coord_15 =  32768;
 constexpr scalar_coordinate_type max_coord_16 =  65536;
 constexpr scalar_coordinate_type max_coord_17 = 131072;
 constexpr scalar_coordinate_type max_coord_18 = 262144;
 constexpr scalar_coordinate_type max_coord_19 = 524288;

 constexpr std::array<scalar_coordinate_type,20> max_coord_arr = {{ max_coord_0,  max_coord_1,  max_coord_2,
                                                                   max_coord_3,  max_coord_4,  max_coord_5,
                                                                   max_coord_6,  max_coord_7,  max_coord_8,
                                                                   max_coord_9,  max_coord_10, max_coord_11,
                                                                   max_coord_12, max_coord_13, max_coord_14,
                                                                   max_coord_15, max_coord_16, max_coord_17,
                                                                   max_coord_18, max_coord_19  }};



}

}

#endif 