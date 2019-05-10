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

//! \file      GDCylinder.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      May 9 2019

#include "GDCylinder.h"

using namespace tinyxml2;

//***************************************************************

GDCylinder::GDCylinder(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, XMLElement *element, const int &physicalEntity, std::string fileName) :
GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement *sousElement(element->FirstChildElement("dataCylinder"));
  if (sousElement == NULL) throw ErrorXMLElement("dataCylinder", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //Radius
  error = sousElement->QueryDoubleAttribute("radius", &m_radius);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("radius", fileName, __FILE__, __LINE__);
  //Axe1
  std::string axe(sousElement->Attribute("axe1"));
  Tools::uppercase(axe);
  if (axe == "X"){ m_axe1 = X; }
  else if (axe == "Y"){ m_axe1 = Y; }
  else if (axe == "Z"){ m_axe1 = Z; }
  else { throw ErrorXMLAttribut("axe1", fileName, __FILE__, __LINE__); }
  //Axe2
  axe = sousElement->Attribute("axe2");
  Tools::uppercase(axe);
  if (axe == "X"){ m_axe2 = X; }
  else if (axe == "Y"){ m_axe2 = Y; }
  else if (axe == "Z"){ m_axe2 = Z; }
  else { throw ErrorXMLAttribut("axe2", fileName, __FILE__, __LINE__); }
  //Length
  error = sousElement->QueryDoubleAttribute("length", &m_length);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("length", fileName, __FILE__, __LINE__);
  if (m_axe1 == X) {
    if (m_axe2 == Y) { m_axe3 = Z; }
    else             { m_axe3 = Y; }
  }
  else if (m_axe1 == Y) {
    if (m_axe2 == X) { m_axe3 = Z; }
    else             { m_axe3 = X; }
  }
  else {
    if (m_axe2 == X) { m_axe3 = Y; }
    else             { m_axe3 = X; }
  }
  //Cylinder center
  double x(0.), y(0.), z(0.);
  XMLElement *center(sousElement->FirstChildElement("center"));
  if (center == NULL) throw ErrorXMLElement("center", fileName, __FILE__, __LINE__);
  error = center->QueryDoubleAttribute("x", &x);
  error = center->QueryDoubleAttribute("y", &y);
  error = center->QueryDoubleAttribute("z", &z);
  m_centerPos.setXYZ(x, y, z);
}

//***************************************************************

GDCylinder::~GDCylinder(){}

//***************************************************************

bool GDCylinder::belong(Coord &posElement, const int &lvl) const
{
  double sum(0.);
  std::vector<Axe> axes;
  axes.push_back(m_axe1);
  axes.push_back(m_axe2);

  for (unsigned int i = 0; i < axes.size(); i++)
  {
    switch (axes[i])
    {
    case X:
      sum += std::pow(posElement.getX() - m_centerPos.getX(), 2.); break;
    case Y:
      sum += std::pow(posElement.getY() - m_centerPos.getY(), 2.); break;
    case Z:
      sum += std::pow(posElement.getZ() - m_centerPos.getZ(), 2.); break;
    }
  }

  double relativePosElement(0.);
  switch (axes[3])
  {
  case X:
    relativePosElement = posElement.getX() - m_centerPos.getX(); break;
  case Y:
    relativePosElement = posElement.getY() - m_centerPos.getY(); break;
  case Z:
    relativePosElement = posElement.getZ() - m_centerPos.getZ(); break;
  }
  
  if ((sum <= m_radius*m_radius) && (std::fabs(relativePosElement) <= m_length/2.)) return true;
  return false;
}