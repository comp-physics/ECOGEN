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

//! \file      GDEllipsoid.cpp
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "GDEllipsoid.h"

using namespace tinyxml2;

//***************************************************************

GDEllipsoid::GDEllipsoid(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, XMLElement *element, const int &physicalEntity, std::string fileName) :
  GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement *sousElement(element->FirstChildElement("dataEllipsoid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataEllipsoid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //radius1
  error = sousElement->QueryDoubleAttribute("radius1", &m_radius1);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("radius1", fileName, __FILE__, __LINE__);
  //radius2
  error = sousElement->QueryDoubleAttribute("radius2", &m_radius2);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("radius2", fileName, __FILE__, __LINE__);
  //radius3
  error = sousElement->QueryDoubleAttribute("radius3", &m_radius3);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("radius3", fileName, __FILE__, __LINE__);
  //Axe1
  std::string axe(sousElement->Attribute("axe1"));
  Tools::uppercase(axe);
  if (axe == "X") { m_axe1 = X; }
  else if (axe == "Y") { m_axe1 = Y; }
  else if (axe == "Z") { m_axe1 = Z; }
  else { throw ErrorXMLAttribut("axe1", fileName, __FILE__, __LINE__); }
  //Axe2
  axe = sousElement->Attribute("axe2");
  Tools::uppercase(axe);
  if (axe == "X") { m_axe2 = X; }
  else if (axe == "Y") { m_axe2 = Y; }
  else if (axe == "Z") { m_axe2 = Z; }
  else { throw ErrorXMLAttribut("axe2", fileName, __FILE__, __LINE__); }
  //Axe3
  axe = sousElement->Attribute("axe3");
  Tools::uppercase(axe);
  if (axe == "X") { m_axe3 = X; }
  else if (axe == "Y") { m_axe3 = Y; }
  else if (axe == "Z") { m_axe3 = Z; }
  else { throw ErrorXMLAttribut("axe3", fileName, __FILE__, __LINE__); }
  //Ellipsoid center
  double x(0.), y(0.), z(0.);
  XMLElement *center(sousElement->FirstChildElement("center"));
  if (center == NULL) throw ErrorXMLElement("center", fileName, __FILE__, __LINE__);
  error = center->QueryDoubleAttribute("x", &x);
  error = center->QueryDoubleAttribute("y", &y);
  error = center->QueryDoubleAttribute("z", &z);
  m_centerPos.setXYZ(x, y, z);
}

//***************************************************************

GDEllipsoid::~GDEllipsoid() {}

//***************************************************************

bool GDEllipsoid::belong(Coord &posElement, const int &lvl) const
{
  double sum(0.);
  std::vector<Axe> axes;
  axes.push_back(m_axe1);
  axes.push_back(m_axe2);
  axes.push_back(m_axe3);

  for (unsigned int i = 0; i < axes.size(); i++)
  {
    switch (axes[i])
    {
    case X:
      if (i == 0)      { sum += std::pow((posElement.getX() - m_centerPos.getX()) / m_radius1, 2.); break; }
      else if (i == 1) { sum += std::pow((posElement.getX() - m_centerPos.getX()) / m_radius2, 2.); break; }
      else             { sum += std::pow((posElement.getX() - m_centerPos.getX()) / m_radius3, 2.); break; }
    case Y:
      if (i == 0)      { sum += std::pow((posElement.getY() - m_centerPos.getY()) / m_radius1, 2.); break; }
      else if (i == 1) { sum += std::pow((posElement.getY() - m_centerPos.getY()) / m_radius2, 2.); break; }
      else             { sum += std::pow((posElement.getY() - m_centerPos.getY()) / m_radius3, 2.); break; }
    case Z:
      if (i == 0)      { sum += std::pow((posElement.getZ() - m_centerPos.getZ()) / m_radius1, 2.); break; }
      else if (i == 1) { sum += std::pow((posElement.getZ() - m_centerPos.getZ()) / m_radius2, 2.); break; }
      else             { sum += std::pow((posElement.getZ() - m_centerPos.getZ()) / m_radius3, 2.); break; }
    }
  }

  if (sum <= 1.) return true;
  return false;
}