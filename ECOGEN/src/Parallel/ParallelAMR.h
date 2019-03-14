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

#ifndef PARALLEL_AMR_H
#define PARALLEL_AMR_H

//! \file      ParallelAMR.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      February 13 2019

#include <mpi.h>
#include "../Tools.h"
#include "../Models/Phase.h"
#include "../Order1/Cell.h"
//KS//BD// Not used for now, see if we can conserve it or erase it
class ParallelAMR : public Parallel 
{
public:
  ParallelAMR();
  virtual ~ParallelAMR();

  //void setNeighbour(const int neighbour);
  //void addElementToSend(const int neighbour, int* numberElement, const int &numberElements);
  //void addElementToReceive(const int neighbour, int* numberElement, const int &numberElements);

};

#endif // PARALLEL_H
