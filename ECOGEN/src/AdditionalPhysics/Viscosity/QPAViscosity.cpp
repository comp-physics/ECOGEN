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

//! \file      QAPViscosity.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "QAPViscosity.h"
#include <iostream>

//***********************************************************************

QAPViscosity::QAPViscosity(){}

//***********************************************************************

QAPViscosity::QAPViscosity(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_grads(3)
{
  m_variableNames.resize(3);
  m_numPhases.resize(3);
  for (int i = 0; i < 3; ++i) {
    m_grads[i] = 0.;
    m_numPhases[i] = -1;
  }
  m_variableNames[0] = "u";
  m_variableNames[1] = "v";
  m_variableNames[2] = "w";
}

//***********************************************************************

QAPViscosity::~QAPViscosity(){}

//***********************************************************************

void QAPViscosity::computeQuantities(Cell* cell)
{
  cell->computeGradient(m_grads, m_variableNames, m_numPhases);
}

//***********************************************************************

void QAPViscosity::setGrad(const Coord &grad, int num)
{
  m_grads[num-1] = grad;
}

//***********************************************************************

Coord QAPViscosity::getGrad(int num) const
{
  return m_grads[num-1];
}

//***********************************************************************