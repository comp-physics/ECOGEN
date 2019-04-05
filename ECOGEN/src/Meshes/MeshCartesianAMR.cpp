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

//! \file      MeshCartesianAMR.cpp
//! \author    K. Schmidmayer, F. Petitpas, B. Dorschner
//! \version   1.0
//! \date      February 19 2019

#include <algorithm>

#include "MeshCartesianAMR.h"

//***********************************************************************

MeshCartesianAMR::MeshCartesianAMR(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
  std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ,
	int lvlMax, double criteriaVar, bool varRho, bool varP, bool varU, bool varAlpha, double xiSplit, double xiJoin) :
  MeshCartesian(lX, numberCellsX, lY, numberCellsY, lZ, numberCellsZ, stretchX, stretchY, stretchZ),
  m_lvlMax(lvlMax), m_criteriaVar(criteriaVar), m_varRho(varRho), m_varP(varP), m_varU(varU), m_varAlpha(varAlpha), m_xiSplit(xiSplit), m_xiJoin(xiJoin)
{
  m_type = AMR;
}

//***********************************************************************

MeshCartesianAMR::~MeshCartesianAMR(){}

//***********************************************************************

int MeshCartesianAMR::initializeGeometrie(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, bool pretraitementParallele, std::string ordreCalcul)
{
  this->meshStretching();
  this->initializeGeometrieAMR(cells, cellsGhost, cellInterfaces, ordreCalcul);
  return m_geometrie;
}


//***********************************************************************

void MeshCartesianAMR::initializeGeometrieAMR(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string ordreCalcul)
{
  int ix, iy, iz;

  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;
  
  //Domain decomposition
  //--------------------
  m_decomp = decomposition::Decomposition({{m_numberCellsXGlobal,m_numberCellsYGlobal,m_numberCellsZGlobal}}, Ncpu);
  auto keys = m_decomp.get_keys(rankCpu);

  for(unsigned int i = 0; i < keys.size(); ++i)
  {    
    if (ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
    else { cells.push_back(new CellO2); }
    m_elements.push_back(new ElementCartesian());
    m_elements[i]->setKey(keys[i]);
    cells[i]->setElement(m_elements[i], i);
  }

  //Create cells and elements
  //-------------------------
  this->assignElementProperties(cells, keys);

  //Create cell interfaces, faces and ghost cells
  //---------------------------------------------
  m_numberCellsCalcul = cells.size();
  createCellInterfacesFacesAndGhostCells(cells, cellsGhost, cellInterfaces, ordreCalcul);
  m_numberCellsTotal = cells.size() + cellsGhost.size();
  m_numberFacesTotal = cellInterfaces.size();

  std::cout<<"cpu "<<rankCpu
    << " m_numberCellsCalcul "<<m_numberCellsCalcul
    << " m_numberCellsTotal "<<m_numberCellsTotal
    << " m_numberFacesTotal "<<m_numberFacesTotal
    <<std::endl; //KS//BD// Erase 3D faces for 1D and 2D simulations
}

//***********************************************************************

void MeshCartesianAMR::assignElementProperties(TypeMeshContainer<Cell *> &cells, std::vector<decomposition::Key<3>> &keys)
{
  int ix, iy, iz;
  double volume(0.);
  for(unsigned int i = 0; i < keys.size(); ++i)
  {
    auto coord = keys[i].coordinate();
    ix = coord.x(); iy = coord.y(); iz = coord.z();
    volume = m_dXi[ix] * m_dYj[iy] * m_dZk[iz];
    cells[i]->getElement()->setVolume(volume);

    //CFL lenght
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[ix]); }
    if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[iy]); }
    if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[iz]); }
    if (m_geometrie > 1) lCFL *= 0.6;

    cells[i]->getElement()->setLCFL(lCFL);
    cells[i]->getElement()->setPos(m_posXi[ix], m_posYj[iy], m_posZk[iz]);
    cells[i]->getElement()->setSize(m_dXi[ix], m_dYj[iy], m_dZk[iz]);
  }
}

//***********************************************************************

void MeshCartesianAMR::createCellInterfacesFacesAndGhostCells(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost,     
TypeMeshContainer<CellInterface*> &cellInterfaces, std::string ordreCalcul)
{
   using coordinate_type = decomposition::Key<3>::coordinate_type;
   std::array<decomposition::Key<3>::coordinate_type,6> offsets;
   std::fill(offsets.begin(), offsets.end(), coordinate_type(0));

   double posX=0, posY=0., posZ=0.;

   for(int d = 0; d < 3; d++)
   {
       offsets[2*d][d] =-1;
       offsets[2*d+1][d] =+1;
   }

   for(unsigned int i = 0; i < cells.size(); ++i)
   {
       const auto coord = cells[i]->getElement()->getKey().coordinate();
       const auto ix = coord.x(), iy = coord.y(), iz = coord.z();
       for(int idx = 0; idx < offsets.size(); idx++)
       {
           const auto offset=offsets[idx];

           posX = m_posXi[ix] + 0.5*m_dXi[ix]*offset[0];
           posY = m_posYj[iy] + 0.5*m_dYj[iy]*offset[1];
           posZ = m_posZk[iz] + 0.5*m_dZk[iz]*offset[2];

           Coord normal, tangent,binormal;
           normal.setXYZ(static_cast<double>(offset[0]), 
                   static_cast<double>(offset[1]), 
                   static_cast<double>(offset[2])); 

           //Xdir
           if (offset[0] == 1) 
           {
               tangent.setXYZ( 0.,1.,0.); 
               binormal.setXYZ(0.,0.,1.); 
           }
           if (offset[0] == -1)
           {
               tangent.setXYZ( 0.,-1.,0.); 
               binormal.setXYZ(0.,0.,1.); 
           }

           //Ydir
           if (offset[1] == 1)
           {
               tangent.setXYZ( -1.,0.,0.); 
               binormal.setXYZ(0.,0.,1.); 
           }
           if (offset[1] == -1)
           {
               tangent.setXYZ( 1.,0.,0.); 
               binormal.setXYZ(0.,0.,1.); 
           }

           //Zdir
           if (offset[2] == 1) 
           {
               tangent.setXYZ( 1.,0.,0.); 
               binormal.setXYZ(0.,1.,0.); 
           }
           if (offset[2] == -1)
           {
               tangent.setXYZ(-1.,0.,0.); 
               binormal.setXYZ(0.,1.,0.); 
           }

           auto neighborCell = cells[i]->getElement()->getKey().coordinate() + offset;
           if (!m_decomp.is_inside(neighborCell)) //Offset is at a physical boundary
           {
               //Create boundary cell interface
               if (offset[0] == 1) //xDir=N
                   m_limXp->creeLimite(cellInterfaces);
               if (offset[0] == -1) //xDir=0
                   m_limXm->creeLimite(cellInterfaces);
               if (offset[1] == 1) //yDir=N
                   m_limYp->creeLimite(cellInterfaces);
               if (offset[1] == -1) //yDir=0
                   m_limYm->creeLimite(cellInterfaces);
               if (offset[2] == 1) //zDir=N
                   m_limZp->creeLimite(cellInterfaces);
               if (offset[2] == -1) //zDir=0
                   m_limZm->creeLimite(cellInterfaces);

               cellInterfaces.back()->initialize(cells[i], nullptr);

               cells[i]->addCellInterface(cellInterfaces.back());
               m_faces.push_back(new FaceCartesian());
               cellInterfaces.back()->setFace(m_faces.back());

               if (offset[0])
               {
                   m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
                   m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
               }
               if (offset[1])
               {
                   m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
                   m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
               }
               if (offset[2])
               {
                   m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
                   m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
               }
               m_faces.back()->setPos(posX, posY, posZ);

           }
           else //Offset is an internal cell (ghost or not)
           {
               //Get neighbor key
               auto nKey = cells[i]->getElement()->getKey().neighbor(offset);
               int neighbour = m_decomp.get_rank(nKey);

               if (offset[0]>0 || offset[1]>0 || offset[2]>0) //Positive offset
               {
                   //Create cell interface
                   if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
                   else { cellInterfaces.push_back(new CellInterfaceO2); }     

                   m_faces.push_back(new FaceCartesian());
                   cellInterfaces.back()->setFace(m_faces.back());

                   if (offset[0])
                   {
                       m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
                       m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
                   }
                   if (offset[1])
                   {
                       m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
                       m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
                   }
                   if (offset[2])
                   {
                       m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
                       m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
                   }
                   m_faces.back()->setPos(posX, posY, posZ);

                   //Try to find the neighbor cell into the non-ghost cells
                   auto it = std::find_if(cells.begin(), cells.end(),
                           [&nKey](Cell* _k0){ 
                           return _k0->getElement()->getKey() == nKey;
                            });

                   if (it != cells.end()) //Neighbor cell is a non-ghost cell
                   {
                       //Update cell interface
                       cellInterfaces.back()->initialize(cells[i], *it);
                       cells[i]->addCellInterface(cellInterfaces.back());
                       (*it)->addCellInterface(cellInterfaces.back());
                   }
                   else //Neighbor cell is a ghost cell
                   {
                       //Try to find the neighbor cell into the already created ghost cells
                       auto it2 = std::find_if(cellsGhost.begin(), cellsGhost.end(),
                         [&nKey](Cell* _k0){ 
                         return _k0->getElement()->getKey() == nKey;
                          });

                       if (it2 == cellsGhost.end()) //Ghost cell does not exist
                       {
                          //Create ghost cell and update cell interface
                         if (ordreCalcul == "FIRSTORDER") { cellsGhost.push_back(new CellGhost); }
                         else { cellsGhost.push_back(new CellO2Ghost); }
                         m_elements.push_back(new ElementCartesian());
                         m_elements.back()->setKey(nKey);
                         cellsGhost.back()->setElement(m_elements.back(), cellsGhost.size()-1);
                         cellsGhost.back()->pushBackSlope();
                         parallel.addSlopesToSend(neighbour);
                         parallel.addSlopesToReceive(neighbour);

                         //Update parallel communications
                         parallel.setNeighbour(neighbour);
                         //Try to find the current cell into the already added cells to send
                         auto cKey = cells[i]->getElement()->getKey();
                         auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                           [&cKey](Cell* _k0){ 
                           return _k0->getElement()->getKey() == cKey;
                            });
                         if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
                         {
                           parallel.addElementToSend(neighbour, cells[i]);
                         }
                         parallel.addElementToReceive(neighbour, cellsGhost.back());
                         cellsGhost.back()->setRankOfNeighborCPU(neighbour);

                         const auto coord = nKey.coordinate();
                         const auto nix = coord.x(), niy = coord.y(), niz = coord.z();

                         const double volume = m_dXi[nix] * m_dYj[niy] * m_dZk[niz];
                         cellsGhost.back()->getElement()->setVolume(volume);

                         double lCFL(1.e10);
                         if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[nix]); }
                         if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[niy]); }
                         if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[niz]); }
                         if (m_geometrie > 1) lCFL *= 0.6;

                         cellsGhost.back()->getElement()->setLCFL(lCFL);
                         cellsGhost.back()->getElement()->setPos(m_posXi[nix], m_posYj[niy], m_posZk[niz]);
                         cellsGhost.back()->getElement()->setSize(m_dXi[nix], m_dYj[niy], m_dZk[niz]);

                         //Update pointers cells <-> cell interfaces
                         cellInterfaces.back()->initialize(cells[i], cellsGhost.back());
                         cells[i]->addCellInterface(cellInterfaces.back());
                         cellsGhost.back()->addCellInterface(cellInterfaces.back());
                       }
                       else { //Ghost cell exists
                         //Update parallel communications
                         //Try to find the current cell into the already added cells to send
                         auto cKey = cells[i]->getElement()->getKey();
                         auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                           [&cKey](Cell* _k0){ 
                           return _k0->getElement()->getKey() == cKey;
                            });
                         if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
                         {
                           parallel.addElementToSend(neighbour, cells[i]);
                         }

                         //Update pointers cells <-> cell interfaces
                         cellInterfaces.back()->initialize(cells[i], *it2);
                         cells[i]->addCellInterface(cellInterfaces.back());
                         (*it2)->addCellInterface(cellInterfaces.back());
                         (*it2)->pushBackSlope();
                         parallel.addSlopesToSend(neighbour);
                         parallel.addSlopesToReceive(neighbour);
                       }
                   }
               }
               else //Negative offset
               {
                   //Try to find the neighbor cell into the non-ghost cells
                   auto it = std::find_if(cells.begin(), cells.end(),
                                             [&nKey](Cell* _k0){ 
                                             return _k0->getElement()->getKey() == nKey;
                                              });

                   if (it == cells.end()) //Neighbor cell is a ghost cell
                   {
                     //Create cell interface related to the ghost cell
                     if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
                     else { cellInterfaces.push_back(new CellInterfaceO2); }     

                     m_faces.push_back(new FaceCartesian());
                     cellInterfaces.back()->setFace(m_faces.back());

                     if (offset[0])
                     {
                         normal.setXYZ( 1.,0.,0.); 
                         tangent.setXYZ( 0.,1.,0.); 
                         binormal.setXYZ(0.,0.,1.); 
                         m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
                         m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
                     }
                     if (offset[1])
                     {
                         normal.setXYZ( 0.,1.,0.); 
                         tangent.setXYZ( -1.,0.,0.); 
                         binormal.setXYZ(0.,0.,1.); 
                         m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
                         m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
                     }
                     if (offset[2])
                     {
                         normal.setXYZ( 0.,0.,1.); 
                         tangent.setXYZ( 1.,0.,0.); 
                         binormal.setXYZ(0.,1.,0.); 
                         m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
                         m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
                     }
                     m_faces.back()->setPos(posX, posY, posZ);

                     //Try to find the neighbor cell into the already created ghost cells
                     auto it2 = std::find_if(cellsGhost.begin(), cellsGhost.end(),
                                            [&nKey](Cell* _k0){ 
                                            return _k0->getElement()->getKey() == nKey;
                                             });

                     if (it2 == cellsGhost.end()) //Ghost cell does not exist
                     {
                       //Create ghost cell
                       if (ordreCalcul == "FIRSTORDER") { cellsGhost.push_back(new CellGhost); }
                       else { cellsGhost.push_back(new CellO2Ghost); }
                       m_elements.push_back(new ElementCartesian());
                       m_elements.back()->setKey(nKey);
                       cellsGhost.back()->setElement(m_elements.back(), cellsGhost.size()-1);
                       cellsGhost.back()->pushBackSlope();
                       parallel.addSlopesToSend(neighbour);
                       parallel.addSlopesToReceive(neighbour);

                       //Update parallel communications
                       parallel.setNeighbour(neighbour);
                       //Try to find the current cell into the already added cells to send
                       auto cKey = cells[i]->getElement()->getKey();
                       auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                         [&cKey](Cell* _k0){ 
                         return _k0->getElement()->getKey() == cKey;
                          });
                       if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
                       {
                         parallel.addElementToSend(neighbour, cells[i]);
                       }
                       parallel.addElementToReceive(neighbour, cellsGhost.back());
                       cellsGhost.back()->setRankOfNeighborCPU(neighbour);

                       const auto coord = nKey.coordinate();
                       const auto nix = coord.x(), niy = coord.y(), niz = coord.z();

                       const double volume = m_dXi[nix] * m_dYj[niy] * m_dZk[niz];
                       cellsGhost.back()->getElement()->setVolume(volume);

                       double lCFL(1.e10);
                       if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[nix]); }
                       if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[niy]); }
                       if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[niz]); }
                       if (m_geometrie > 1) lCFL *= 0.6;

                       cellsGhost.back()->getElement()->setLCFL(lCFL);
                       cellsGhost.back()->getElement()->setPos(m_posXi[nix], m_posYj[niy], m_posZk[niz]);
                       cellsGhost.back()->getElement()->setSize(m_dXi[nix], m_dYj[niy], m_dZk[niz]);

                       //Update pointers cells <-> cell interfaces
                       cellInterfaces.back()->initialize(cellsGhost.back(), cells[i]);
                       cells[i]->addCellInterface(cellInterfaces.back());
                       cellsGhost.back()->addCellInterface(cellInterfaces.back());
                     }
                     else //Ghost cell exists
                     {
                       //Update parallel communications
                       //Try to find the current cell into the already added cells to send
                       auto cKey = cells[i]->getElement()->getKey();
                       auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                         [&cKey](Cell* _k0){ 
                         return _k0->getElement()->getKey() == cKey;
                          });
                       if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
                       {
                         parallel.addElementToSend(neighbour, cells[i]);
                       }

                       //Update pointers cells <-> cell interfaces
                       cellInterfaces.back()->initialize(*it2, cells[i]);
                       cells[i]->addCellInterface(cellInterfaces.back());
                       (*it2)->addCellInterface(cellInterfaces.back());
                       (*it2)->pushBackSlope();
                       parallel.addSlopesToSend(neighbour);
                       parallel.addSlopesToReceive(neighbour);
                     }
                   }
               } //Negative offset
           } //Offset is an internal cell (ghost or not)
       } //Offsets
   } //Internal, non-ghost cells

   if (Ncpu > 1) {
     for(int i=0;i<Ncpu;++i)
     {
      std::sort(parallel.getElementsToReceive(i).begin(),parallel.getElementsToReceive(i).end(),[&]( Cell* child0, Cell* child1 )
      {
        return child0->getElement()->getKey()< child1->getElement()->getKey();
      });
     }
   }
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinementInitialization(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, TypeMeshContainer<CellInterface *> *cellInterfacesLvl,
  const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, std::vector<GeometricalDomain*> &domains,
  Eos **eos, const int &resumeSimulation, std::string ordreCalcul, const int &numberPhases, const int &numberTransports)
{
  nbCellsTotalAMR = m_numberCellsCalcul;

  if (resumeSimulation == 0) { //Only for simulation from input files
    for (int iterInit = 0; iterInit < 2; iterInit++) {
      for (int lvl = 0; lvl < m_lvlMax; lvl++) {
        if (Ncpu > 1) { parallel.communicationsPrimitives(eos, lvl); }
        this->procedureRaffinement(cellsLvl, cellsLvlGhost, cellInterfacesLvl, lvl, addPhys, model, nbCellsTotalAMR, eos);
        //if (Ncpu > 1) { this->parallelLoadBalancingAMR(cellsLvl, cellsLvlGhost, cellInterfacesLvl, ordreCalcul, numberPhases, numberTransports, addPhys, model); } //KS//BD//
        this->parallelLoadBalancingAMR(cellsLvl, cellsLvlGhost, cellInterfacesLvl, ordreCalcul, numberPhases, numberTransports, addPhys, model);
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->fill(domains, m_lvlMax);
        }
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->completeFulfillState();
        }
        for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
          cellsLvl[lvl][i]->averageChildrenInParent();
        }
      }
    }
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      if (Ncpu > 1) { parallel.communicationsPrimitives(eos, lvl); }
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->completeFulfillState(); }
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinement(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, TypeMeshContainer<CellInterface *> *cellInterfacesLvl, const int &lvl,
  const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, Eos **eos)
{
  //1) Calcul de Xi dans chaque cell de niveau lvl
  //-------------------------------------------------
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroXi(); }
  for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, m_varP, m_varU, m_varAlpha); }
  //bool varP2(true);
  //if (lvl >= 5) { varP2 = false; }
  //for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, varP2, m_varU, m_varAlpha); }
  //for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
  //  double x(0.), y(0.), z(0.);
  //  x = cellsLvl[lvl][i]->getPosition().getX();
  //  y = cellsLvl[lvl][i]->getPosition().getY();
  //  //z = cellsLvl[lvl][i]->getPosition().getZ();
  //  //if (std::pow((x*x + y*y + z*z), 0.5) > 500.e-6) {
  //  //if (std::pow((x*x + y*y), 0.5) > 6.e-4) {
  //  //if ((x > 250e-6) || (y > 200.e-6)) {
  //  //if (x > 15.) {
  //  if (std::pow((x*x + y * y), 0.5) > 5.) {
  //      cellsLvl[lvl][i]->setToZeroXi();
  //  }
  //}
  if (Ncpu > 1) { parallel.communicationsXi( lvl); }
  
  //2) Smoothing de Xi
  //------------------
  for (int iterDiff = 0; iterDiff < 2; iterDiff++) { //Arbitrary number of iterations
		//Mise a zero cons xi
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroConsXi(); }

    //Calcul des "flux"
    for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeFluxXi(); }

    //Evolution temporelle
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->timeEvolutionXi(); }
		if (Ncpu > 1) { parallel.communicationsXi( lvl); }
  }

	if (lvl < m_lvlMax) {
    int lvlPlus1 = lvl + 1;
    //3) Raffinement des cells et cell interfaces
    //-------------------------------------------
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseRefine(m_xiSplit, m_numberCellsY, m_numberCellsZ, addPhys, model, nbCellsTotalAMR); }

    //4) Deraffinement des cells et cell interfaces
    //---------------------------------------------
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseUnrefine(m_xiJoin, nbCellsTotalAMR); }

    if (Ncpu > 1) {
      //5) Raffinement et deraffinement des cells fantomes
      //-----------------------------------------------------
      //Communication split + Raffinement et deraffinement des cells fantomes + Reconstruction du tableau de cells fantomes de niveau lvl + 1
      parallel.communicationsSplit(lvl);
      cellsLvlGhost[lvlPlus1].clear();
      for (unsigned int i = 0; i < cellsLvlGhost[lvl].size(); i++) { cellsLvlGhost[lvl][i]->chooseRefineDeraffineGhost(m_numberCellsY, m_numberCellsZ, addPhys, model, cellsLvlGhost); }
      //Communications primitives pour mettre a jour les cells deraffinees
      parallel.communicationsPrimitives(eos, lvl);

      //6) Mise a jour des communications persistantes au niveau lvl + 1
      //----------------------------------------------------------------
      parallel.communicationsNumberGhostCells(lvlPlus1);	//Communication des numbers d'elements a envoyer et a recevoir de chaque cote de la limite parallele
      parallel.updatePersistentCommunicationsLvl(lvlPlus1, m_geometrie);
    }

    //7) Reconstruction des tableaux de cells et cell interfaces lvl + 1
    //------------------------------------------------------------------
    cellsLvl[lvlPlus1].clear();
    cellInterfacesLvl[lvlPlus1].clear();
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalCellInterfacesArrays(cellsLvl, cellInterfacesLvl); }
    for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->constructionTableauCellInterfacesExternesLvl(cellInterfacesLvl); }
  }
}

//***********************************************************************

std::string MeshCartesianAMR::whoAmI() const
{
  return "CARTESIAN_AMR";
}

//**************************************************************************
//******************************** PRINTING ********************************
//**************************************************************************

void MeshCartesianAMR::ecritHeaderPiece(std::ofstream &fileStream, TypeMeshContainer<Cell *> *cellsLvl) const
{
  int numberCells = 0, numberPointsParMaille = 4;
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        numberCells += 1;
      }
    }
  }
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }

  fileStream << "    <Piece NumberOfPoints=\"" << numberPointsParMaille*numberCells << "\" NumberOfCells=\"" << numberCells << "\">" << std::endl;
}

//***********************************************************************

void MeshCartesianAMR::recupereNoeuds(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int dimZ = 0;
  if (m_numberCellsZ > 1) dimZ = 1;

  double dXsur2(0.), dYsur2(0.), dZsur2(0.);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    dXsur2 = 0.; dYsur2 = 0.; dZsur2 = 0.;
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        dXsur2 = 0.5*cellsLvl[lvl][i]->getSizeX();
        dYsur2 = 0.5*cellsLvl[lvl][i]->getSizeY();
        dZsur2 = 0.5*cellsLvl[lvl][i]->getSizeZ();
        //Point 0
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 1
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 2
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 3
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);

        if (dimZ > 0.99) {
          //Point 4
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 5
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 6
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 7
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
        }
      } //Fin cell non split
    } //Fin Cells
  } //Fin Levels
}

//***********************************************************************

void MeshCartesianAMR::recupereConnectivite(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int dimZ(0);
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { dimZ = 1; numberPointsParMaille = 8; }

  if (dimZ < 0.99) {
    int numCell(0);
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) {
          jeuDonnees.push_back(numCell*numberPointsParMaille);
          jeuDonnees.push_back(numCell*numberPointsParMaille+1);
          jeuDonnees.push_back(numCell*numberPointsParMaille+2);
          jeuDonnees.push_back(numCell*numberPointsParMaille+3);
          numCell++;
        }
      }
    }
  }
  else {
    int numCell(0);
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) {
          jeuDonnees.push_back(numCell*numberPointsParMaille);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 1);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 2);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 3);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 4);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 5);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 6);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 7);
          numCell++;
        }
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereOffsets(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }
  int numCell(0);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        jeuDonnees.push_back((numCell + 1)*numberPointsParMaille);
        numCell++;
      }
    }
  }
}

//****************************************************************************

void MeshCartesianAMR::recupereTypeCell(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int type(9);
  if (m_numberCellsZ > 1) { type = 12; }
  int numCell(0);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        jeuDonnees.push_back(type);
        numCell++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereDonnees(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase) const
{
  jeuDonnees.clear();
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        if (var > 0) { //On veut recuperer les donnees scalars
          if (phase >= 0) { jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnScalar(var)); }      //Donnees de phases
          else if (phase == -1) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnScalar(var)); }   //Donnees de mixture
          else if (phase == -2) { jeuDonnees.push_back(cellsLvl[lvl][i]->getTransport(var - 1).getValue()); }
          else if (phase == -3) { jeuDonnees.push_back(cellsLvl[lvl][i]->getXi()); }
          else if (phase == -4) { jeuDonnees.push_back(cellsLvl[lvl][i]->getGradient()); }
          else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Donnees de phases
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getZ());
          }
          else if(phase == -1){  //Donnees de mixture
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getZ());
          }
          else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
        } //Fin vecteur
      } //Fin split
    } //fin lvl
  } //fin levels
}

//****************************************************************************

void MeshCartesianAMR::setDataSet(std::vector<double> &jeuDonnees, TypeMeshContainer<Cell *> *cellsLvl, const int var, int phase) const
{
  int iterDataSet(0);
  Coord vec;
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        if (var > 0) { //Scalars data are first set
          if (phase >= 0) { cellsLvl[lvl][i]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
          else if (phase == -1) { cellsLvl[lvl][i]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
          else if (phase == -2) { cellsLvl[lvl][i]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
          else if (phase == -3) { cellsLvl[lvl][i]->setXi(jeuDonnees[iterDataSet++]); } //xi indicator
          else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Phases data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[lvl][i]->getPhase(phase)->setVector(-var, vec);
            iterDataSet += 3;
          }
          else if (phase == -1) {  //Mixture data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[lvl][i]->getMixture()->setVector(-var, vec);
            iterDataSet += 3;
          }
          else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
        } //Fin vecteur
      } // Fin split
    } // Fin lvl
  } // Fin levels
}

//***********************************************************************

void MeshCartesianAMR::refineCell(Cell *cell, const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR)
{
  cell->refineCellAndCellInterfaces(m_numberCellsY, m_numberCellsZ, addPhys, model);
  nbCellsTotalAMR += cell->getNumberCellsChildren() - 1;
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void MeshCartesianAMR::initializePersistentCommunications(const int numberPhases, const int numberTransports, const TypeMeshContainer<Cell *> &cells, std::string ordreCalcul)
{
	m_numberPhases = numberPhases;
	m_numberTransports = numberTransports;
	int numberVariablesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedVariables();
	numberVariablesPhaseATransmettre *= m_numberPhases;
	int numberVariablesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedVariables();
	int m_numberPrimitiveVariables = numberVariablesPhaseATransmettre + numberVariablesMixtureATransmettre + m_numberTransports;
  int m_numberSlopeVariables(0);
  if (ordreCalcul == "SECONDORDER") {
    int numberSlopesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedSlopes();
    numberSlopesPhaseATransmettre *= m_numberPhases;
    int numberSlopesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedSlopes();
    m_numberSlopeVariables = numberSlopesPhaseATransmettre + numberSlopesMixtureATransmettre + m_numberTransports + 1 + 1; //+1 for the interface detection + 1 for slope index
  }
	parallel.initializePersistentCommunicationsAMR(m_numberPrimitiveVariables, m_numberSlopeVariables, m_numberTransports, m_geometrie, m_lvlMax);
}

//***********************************************************************

void MeshCartesianAMR::finalizeParallele(const int &lvlMax)
{
	parallel.finalizeAMR(lvlMax);
}

//***********************************************************************

void MeshCartesianAMR::parallelLoadBalancingAMR(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, TypeMeshContainer<CellInterface *> *cellInterfacesLvl, std::string ordreCalcul,
  const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model)
{
  int iter(0);
  double idealLoadEndPosition(0.);
  double *loadPerCPU = new double[Ncpu];
  typename decomposition::Key<3>::value_type *startPerCPU = new typename decomposition::Key<3>::value_type[Ncpu];
  MPI_Request req_sendNeighborP1;
  MPI_Request req_recvNeighborM1;
  MPI_Request req_sendNeighborP1_2;
  MPI_Request req_recvNeighborM1_2;
  MPI_Status status;

  //LOOP over the following until global balance
  //while (CONDITION)

  //Typical example of one situation for CPU 1 (during 1 balance iteration):
  //
  //                CPU 0       CPU 1                                 CPU 2                           
  // Current load |-------|---------------|-----------------------------------------------------------|
  //                      |   localLoad   |
  //                      |      localLoadEndPosition
  //             localLoadStartPosition
  //
  //                          CPU 0                       CPU 1                       CPU 2
  // Ideal load   |---------------------------|---------------------------|---------------------------|
  //                                 idealLoadStartPosition      idealLoadEndPosition
  //
  // Ideal load shifts    |------------------->
  //                                      |------------------------------->
  //                       idealLoadShiftStart    idealLoadShiftEnd
  //
  // Possible load shifts |--------------->
  //                                      |------------------------------->
  //                    possibleLoadShiftStart   possibleLoadShiftEnd
  //
  //                        CPU 0                       CPU 1                         CPU 2
  // Final load   |-----------------------|-------------------------------|---------------------------|
  //                                      |        finalLocalLoad         |
  //                                      |                   finalLocalLoadEndPosition
  //                          finalLocalLoadStartPosition
  //
  // Note that the possible load shift start gives birth of a new start that may
  // constrain the possible load shift end (situation not present in this example).

  //1) Compute and communicate loads and positions
  //----------------------------------------------

  // //Compute local load
  // double localLoad(0.);
  // for (unsigned int i = 0; i < cellsLvl[0].size(); i++) {
  //   cellsLvl[0][i]->computeLoad(localLoad);
  // }

  // //Communicate overall loads
  // MPI_Allgather(&localLoad, 1, MPI_DOUBLE, loadPerCPU, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  // //Compute ideal load end position (only for the first iteration)
  // if (iter == 0) {
  //   for (int i = 0; i < Ncpu; ++i) { idealLoadEndPosition += loadPerCPU[i]; }
  //   std::cout<<"cpu "<<rankCpu<<" totalLoad "<<idealLoadEndPosition<<std::endl; //KS//BD//
  //   idealLoadEndPosition *= static_cast<double>(rankCpu + 1) / Ncpu;
  //   std::cout<<"cpu "<<rankCpu<<" idealLoadEndPosition "<<idealLoadEndPosition<<std::endl; //KS//BD//
  // }

  // //Compute local load end position
  // double localLoadEndPosition(0.);
  // for (int i = 0; i <= rankCpu; ++i) { localLoadEndPosition += loadPerCPU[i]; }
  // std::cout<<"cpu "<<rankCpu<<" localLoadEndPosition "<<localLoadEndPosition<<std::endl; //KS//BD//

  // //2) Compute and communicate ideal shifts
  // //---------------------------------------

  // //Determine load I should send/receive to/from/ CPU P1, using end positions of local load and ideal load
  // //if (diff > 0), I wish to receive loads from CPU P1
  // //else,          I wish to send loads to CPU P1
  // double idealLoadShiftEnd(0.);
  // idealLoadShiftEnd = idealLoadEndPosition - localLoadEndPosition;

  // //Communicate what I wish to send/receive to/from neighbours
  // double idealLoadShiftStart(0.);
  // if (rankCpu != Ncpu - 1) MPI_Isend(&idealLoadShiftEnd, 1, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1);
  // if (rankCpu != 0)        MPI_Irecv(&idealLoadShiftStart, 1, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1);

  // if (rankCpu != Ncpu - 1) MPI_Wait(&req_sendNeighborP1, &status);
  // if (rankCpu != 0)        MPI_Wait(&req_recvNeighborM1, &status);

  // //3) Compute and communicate possible shifts (real balance)
  // //---------------------------------------------------------

  // //Determine and communicate what I can send/receive to/from neighbours (limited by current local cells/load)
  // double possibleLoadShiftStart(0.), possibleLoadShiftEnd(0.);
  // int newStart(0);
  // int numberOfCellsToSendStart(0), numberOfCellsToSendEnd(0);
  // int numberOfCellsToReceiveStart(0), numberOfCellsToReceiveEnd(0);

  // //For load shift start
  // if (idealLoadShiftStart > 0.) {
  //   //Determine and send possible load shift start
  //   for (unsigned int i = 0; i < cellsLvl[0].size(); i++) {
  //     if (possibleLoadShiftStart >= idealLoadShiftStart) {
  //       newStart = i;
  //       break;
  //     }
  //     cellsLvl[0][i]->computeLoad(possibleLoadShiftStart);
  //     numberOfCellsToSendStart = i + 1;
  //   }
  //   if (rankCpu != 0)        MPI_Isend(&possibleLoadShiftStart, 1, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1);
  //   if (rankCpu != 0)        MPI_Isend(&numberOfCellsToSendStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1_2);
  // }
  // else {
  //   //Receive possible load shift start
  //   if (rankCpu != 0)        MPI_Irecv(&possibleLoadShiftStart, 1, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1);
  //   if (rankCpu != 0)        MPI_Irecv(&numberOfCellsToReceiveStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1_2);
  // }

  // //For load shift end
  // if (idealLoadShiftEnd < 0.) {
  //   //Determine and send possible load shift end
  //   for (int i = cellsLvl[0].size()-1; i >= newStart; i--) {
  //     if (possibleLoadShiftEnd >= (idealLoadShiftEnd)) {
  //       break;
  //     }
  //     cellsLvl[0][i]->computeLoad(possibleLoadShiftEnd);
  //     numberOfCellsToSendEnd = cellsLvl[0].size() - i;
  //   }
  //   if (rankCpu != Ncpu - 1) MPI_Isend(&possibleLoadShiftEnd, 1, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1);
  //   if (rankCpu != Ncpu - 1) MPI_Isend(&numberOfCellsToSendEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1_2);
  // }
  // else {
  //   //Receive possible load shift end
  //   if (rankCpu != Ncpu - 1) MPI_Irecv(&possibleLoadShiftEnd, 1, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1);
  //   if (rankCpu != Ncpu - 1) MPI_Irecv(&numberOfCellsToReceiveEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1_2);
  // }

  // if (rankCpu != Ncpu - 1) MPI_Wait(&req_sendNeighborP1, &status);
  // if (rankCpu != 0)        MPI_Wait(&req_recvNeighborM1, &status);
  // if (rankCpu != Ncpu - 1) MPI_Wait(&req_sendNeighborP1_2, &status);
  // if (rankCpu != 0)        MPI_Wait(&req_recvNeighborM1_2, &status);

  // //Final load......
  // std::cout<<"cpu "<<rankCpu
  // <<" idealLoadShiftStart "<<idealLoadShiftStart<<" idealLoadShiftEnd "<<idealLoadShiftEnd
  // <<" possibleLoadShiftStart "<<possibleLoadShiftStart<<" possibleLoadShiftEnd "<<possibleLoadShiftEnd
  // <<std::endl; //KS//BD//
  // double localLoadStartPosition = localLoadEndPosition - loadPerCPU[rankCpu];
  // localLoadStartPosition += possibleLoadShiftStart;
  // localLoadEndPosition += possibleLoadShiftEnd;
  // double finalLocalLoad =localLoadEndPosition-localLoadStartPosition;
  // std::cout<<"cpu "<<rankCpu<<" localLoadStartPosition "<<localLoadStartPosition<<" localLoadEndPosition "<<localLoadEndPosition
  // <<" initialLocalLoad "<<loadPerCPU[rankCpu]
  // <<" finalLocalLoad "<<finalLocalLoad<<std::endl; //KS//BD//

  // //4) Send/Receive keys of base cells (lvl = 0)
  // //--------------------------------------------
  // std::cout<<"cpu "<<rankCpu
  // <<" numberOfCellsToSendStart "<<numberOfCellsToSendStart<<" numberOfCellsToSendEnd "<<numberOfCellsToSendEnd
  // <<" numberOfCellsToReceiveStart "<<numberOfCellsToReceiveStart<<" numberOfCellsToReceiveEnd "<<numberOfCellsToReceiveEnd
  // <<std::endl; //KS//BD//

  // //Send keys of the base cells
  // if (numberOfCellsToSendStart > 0) {
  //   typename decomposition::Key<3>::value_type *keys = new typename decomposition::Key<3>::value_type[numberOfCellsToSendStart];
  //   for (unsigned int i = 0; i < numberOfCellsToSendStart; i++) {
  //     keys[i] = cellsLvl[0][i]->getElement()->getKey().getIndex();
  //   }
  //   if (rankCpu != 0)        MPI_Isend(keys, numberOfCellsToSendStart, MPI_UNSIGNED_LONG_LONG, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1);
  //   if (rankCpu != 0)        MPI_Wait(&req_recvNeighborM1, &status);
  //   delete[] keys;
  // }

  // if (numberOfCellsToSendEnd > 0) {
  //   typename decomposition::Key<3>::value_type *keys = new typename decomposition::Key<3>::value_type[numberOfCellsToSendEnd];
  //   for (int i = cellsLvl[0].size()-1; i >= cellsLvl[0].size()-1-numberOfCellsToSendEnd; --i) {
  //     keys[i] = cellsLvl[0][i]->getElement()->getKey().getIndex();
  //   }
  //   if (rankCpu != Ncpu - 1) MPI_Isend(keys, numberOfCellsToSendEnd, MPI_UNSIGNED_LONG_LONG, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1);
  //   if (rankCpu != Ncpu - 1) MPI_Wait(&req_sendNeighborP1, &status);
  //   delete[] keys;
  // }

  // //Receive keys of the base cells
  // typename decomposition::Key<3>::value_type *receivedIndicesStart = new typename decomposition::Key<3>::value_type[numberOfCellsToReceiveStart];
  // if (numberOfCellsToReceiveStart > 0) {
  //   if (rankCpu != 0)        MPI_Irecv(receivedIndicesStart, numberOfCellsToReceiveStart, MPI_UNSIGNED_LONG_LONG, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_recvNeighborM1);
  //   if (rankCpu != 0)        MPI_Wait(&req_recvNeighborM1, &status);
  // }

  // typename decomposition::Key<3>::value_type *receivedIndicesEnd = new typename decomposition::Key<3>::value_type[numberOfCellsToReceiveEnd];
  // if (numberOfCellsToReceiveEnd > 0) {
  //   if (rankCpu != Ncpu - 1) MPI_Irecv(receivedIndicesEnd, numberOfCellsToReceiveEnd, MPI_UNSIGNED_LONG_LONG, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_sendNeighborP1);
  //   if (rankCpu != Ncpu - 1) MPI_Wait(&req_sendNeighborP1, &status);
  // }

  // //5) Create the corresponding base cells (received)
  // //-------------------------------------------------

  // //Create and insert cells and elements
  // TypeMeshContainer<Cell *> bufferCellsStart(numberOfCellsToReceiveStart);
  // TypeMeshContainer<Cell *> bufferCellsEnd(numberOfCellsToReceiveEnd);
  // TypeMeshContainer<Cell *> bufferCellsNew;
  // std::vector<decomposition::Key<3>> receivedKeysStart(numberOfCellsToReceiveStart);
  // std::vector<decomposition::Key<3>> receivedKeysEnd(numberOfCellsToReceiveEnd);
  // for(unsigned int i = 0; i < numberOfCellsToReceiveStart; ++i) {
  //   if (ordreCalcul == "FIRSTORDER") { bufferCellsStart[i] = new Cell; }
  //   else { bufferCellsStart[i] = new CellO2; }
  //   m_elements.push_back(new ElementCartesian());
  //   receivedKeysStart[i] = decomposition::Key<3>(receivedIndicesStart[i]);
  //   m_elements.back()->setKey(receivedKeysStart[i]);
  //   bufferCellsStart[i]->setElement(m_elements.back(), i);
  // }
  // for(unsigned int i = 0; i < numberOfCellsToReceiveEnd; ++i) {
  //   if (ordreCalcul == "FIRSTORDER") { bufferCellsEnd[i] = new Cell; }
  //   else { bufferCellsEnd[i] = new CellO2; }
  //   m_elements.push_back(new ElementCartesian());
  //   receivedKeysEnd[i] = decomposition::Key<3>(receivedIndicesEnd[i]);
  //   m_elements.back()->setKey(receivedKeysEnd[i]);
  //   bufferCellsEnd[i]->setElement(m_elements.back(), i);
  // }
  // delete[] receivedIndicesStart;
  // delete[] receivedIndicesEnd;

  // cellsLvl[0].insert(cellsLvl[0].begin(), bufferCellsStart.begin(), bufferCellsStart.end());
  // cellsLvl[0].insert(cellsLvl[0].end(), bufferCellsEnd.begin(), bufferCellsEnd.end());

  // //Assigning element properties
  // this->assignElementProperties(bufferCellsStart, receivedKeysStart);
  // this->assignElementProperties(bufferCellsEnd, receivedKeysEnd);

  // bufferCellsNew.insert(bufferCellsNew.end(), bufferCellsStart.begin(), bufferCellsStart.end());
  // bufferCellsNew.insert(bufferCellsNew.end(), bufferCellsEnd.begin(), bufferCellsEnd.end());

  // //6) Erase the pointers to the corresponding base cells from cellsLvl (cells are not deleted yet)
  // //-----------------------------------------------------------------------------------------------
  // TypeMeshContainer<Cell *> temporaryCellsStart(numberOfCellsToSendStart);
  // TypeMeshContainer<Cell *> temporaryCellsEnd(numberOfCellsToSendEnd);
  // for(unsigned int i = 0; i < numberOfCellsToSendStart; ++i) {
  //   temporaryCellsStart[i] = cellsLvl[0][i];
  // }
  // for(unsigned int i = 0; i < numberOfCellsToSendEnd; ++i) {
  //   temporaryCellsEnd[i] = cellsLvl[0][cellsLvl[0].size()-1-i];
  // }
  // cellsLvl[0].erase(cellsLvl[0].begin(), cellsLvl[0].begin()+numberOfCellsToSendStart);
  // cellsLvl[0].erase(cellsLvl[0].end()-numberOfCellsToSendEnd, cellsLvl[0].end());

  // //7) Update the domain decomposition with new starts and ends
  // //-----------------------------------------------------------
  // typename decomposition::Key<3>::value_type localStart = cellsLvl[0][0]->getElement()->getKey().getIndex();
  // MPI_Allgather(&localStart, 1, MPI_UNSIGNED_LONG_LONG, startPerCPU, 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
  // for( int i = 0; i < Ncpu; ++i) {
  //   decomposition::Key<3> key_tmp(startPerCPU[i]);
  //   m_decomp.start_key(i)=key_tmp;
  // }

  //8) Create cell interfaces, faces and ghost cells level 0
  //--------------------------------------------------------
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " before balance"
// <<std::endl;
// for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
// std::cout<<"cpu "<<rankCpu
// << " cellsLvl["<<lvl<<"].size() "<<cellsLvl[lvl].size()
// << " cellInterfacesLvl["<<lvl<<"].size() "<<cellInterfacesLvl[lvl].size()
// <<std::endl;
// }
  for (int b = 0; b < cellInterfacesLvl[0].size(); b++) { delete cellInterfacesLvl[0][b]; }
  for (int i = 0; i < cellsLvl[0].size(); i++) { cellsLvl[0][i]->clearExternalCellInterfaces(m_numberCellsY, m_numberCellsZ); }
  for (int i = 0; i < cellsLvlGhost[0].size(); i++) { delete cellsLvlGhost[0][i]; }
  cellInterfacesLvl[0].clear();
  cellsLvlGhost[0].clear();
  //parallel.clearElementsAndSlopesToSendAndReceive(); //KS//BD// Comment for serial test
  m_numberCellsCalcul = cellsLvl[0].size();
  createCellInterfacesFacesAndGhostCells(cellsLvl[0], cellsLvlGhost[0], cellInterfacesLvl[0], ordreCalcul);
  m_numberCellsTotal = cellsLvl[0].size() + cellsLvlGhost[0].size();
  m_numberFacesTotal = cellInterfacesLvl[0].size();

  std::cout<<"cpu "<<rankCpu
    << " m_numberCellsCalcul "<<m_numberCellsCalcul
    << " m_numberCellsTotal "<<m_numberCellsTotal
    << " m_numberFacesTotal "<<m_numberFacesTotal
    <<std::endl;

  //9) Allocate physical variables of cells and cell interfaces level 0
  //-------------------------------------------------------------------
  // for (int i = 0; i < bufferCellsNew.size(); i++) { bufferCellsNew[i]->allocate(numberPhases, numberTransports, addPhys, model); }
  for (int i = 0; i < cellsLvlGhost[0].size(); i++) { cellsLvlGhost[0][i]->allocate(numberPhases, numberTransports, addPhys, model); }
  //Attribution model and slopes to faces
  int allocateSlopeLocal = 1;
  for (int b = 0; b < cellInterfacesLvl[0].size(); b++) {
    cellInterfacesLvl[0][b]->associeModel(model);
    cellInterfacesLvl[0][b]->allocateSlopes(numberPhases, numberTransports, allocateSlopeLocal);
  }

  //Refine external cell interfaces for already splitted cells (the ones not shift around)
  int dim(1);
  if (m_numberCellsZ != 1) { dim = 3; }
  else if (m_numberCellsY != 1) { dim = 2; }
  for (int lvl = 0; lvl < m_lvlMax; lvl++) {
    for (int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (cellsLvl[lvl][i]->getSplit()) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " posX "<<cellsLvl[lvl][i]->getPosition().getX()
// << " posY "<<cellsLvl[lvl][i]->getPosition().getY()
// << " getCellInterfacesSize "<<cellsLvl[lvl][i]->getCellInterfacesSize()
// <<std::endl;
        for (int b = 0; b < cellsLvl[lvl][i]->getCellInterfacesSize(); b++) {
// if (lvl == 1) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " posX "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getX()
// << " posY "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getY()
// <<" lvl "<<cellsLvl[lvl][i]->getCellInterface(b)->getLvl()
// <<std::endl;
// }
          if (!cellsLvl[lvl][i]->getCellInterface(b)->getSplit()) {
// if (lvl == 1) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// <<" getCellInterface "<<b
// <<" whoAmI "<<cellsLvl[lvl][i]->getCellInterface(b)->whoAmI()
// <<std::endl;
// }
// if (lvl == 1) {
// std::cout<<"HERE"<<std::endl;
// }
            cellsLvl[lvl][i]->getCellInterface(b)->raffineCellInterfaceExterne(m_numberCellsY, m_numberCellsZ, 
              cellsLvl[lvl][i]->getElement()->getSizeX(), cellsLvl[lvl][i]->getElement()->getSizeY(), cellsLvl[lvl][i]->getElement()->getSizeZ(), cellsLvl[lvl][i], dim);
          }
// if (cellsLvl[lvl][i]->getCellInterface(b)->getLvl() == 1) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// <<" getCellInterface "<<b
// <<" whoAmI "<<cellsLvl[lvl][i]->getCellInterface(b)->whoAmI()
// <<" parentCell "<<cellsLvl[lvl][i]->getPosition().getX()
// <<" parentCell "<<cellsLvl[lvl][i]->getPosition().getY()
// <<" lvl "<<cellsLvl[lvl][i]->getLvl()
// <<std::endl;
// if (cellsLvl[lvl][i]->getCellInterface(b)->whoAmI() == 0) {
// for (int c = 0; c < cellsLvl[lvl][i]->getCellInterface(b)->getNumberCellInterfacesChildren(); c++) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " posX "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getFace()->getPos().getX()
// << " posY "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getFace()->getPos().getY()
// <<" lvl "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getLvl()
// << " left "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getCellGauche()->getPosition().getX()
// << " "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getCellGauche()->getPosition().getY()
// << " right "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getCellDroite()->getPosition().getX()
// << " "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellInterfaceChild(c)->getCellDroite()->getPosition().getY()
// <<std::endl;
// }
// }
// }
// else if (cellsLvl[lvl][i]->getCellInterface(b)->getLvl() == 2) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " posX "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getX()
// << " posY "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getY()
// <<" lvl "<<cellsLvl[lvl][i]->getCellInterface(b)->getLvl()
// << " left "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellGauche()->getPosition().getX()
// << " "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellGauche()->getPosition().getY()
// << " right "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellDroite()->getPosition().getX()
// << " "<<cellsLvl[lvl][i]->getCellInterface(b)->getCellDroite()->getPosition().getY()
// <<std::endl;
// }
        }
      }
    }
    //Reconstruction of the arrays of cells and cell interfaces of lvl + 1
    if (lvl < m_lvlMax) {
      cellsLvl[lvl + 1].clear();
      cellInterfacesLvl[lvl + 1].clear();
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalCellInterfacesArrays(cellsLvl, cellInterfacesLvl); }
      for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->constructionTableauCellInterfacesExternesLvl(cellInterfacesLvl); }
    }
  }

// int lvl =2; //KS//BD//
// for (int i = 0; i < cellsLvl[lvl].size(); i++) {
// for (int b = 0; b < cellsLvl[lvl][i]->getCellInterfacesSize(); b++) {
// if (lvl == 2) {
// std::cout<<"cpu "<<rankCpu //KS//BD//
// << " posX "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getX()
// << " posY "<<cellsLvl[lvl][i]->getCellInterface(b)->getFace()->getPos().getY()
// <<" lvl "<<cellsLvl[lvl][i]->getCellInterface(b)->getLvl()
// <<std::endl;
// }
// }
// }

  //10) Send/Receive physical values
  //--------------------------------
  // for (int i = 0; i < bufferCellsNew.size(); i++) { bufferCellsNew[i]->fill(domains, m_lvlMax); }
  // for (int i = 0; i < cellsLvlGhost[0].size(); i++) { cellsLvlGhost[0][i]->fill(domains, m_lvlMax); }



  // //Complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  // for (int i = 0; i < bufferCellsNew.size(); i++) { bufferCellsNew[i]->completeFulfillState(); }



  // //7) Intialization of persistant communications for parallel computing
  // //--------------------------------------------------------------------
  // m_mesh->initializePersistentCommunications(m_numberPhases, m_numberTransports, m_cellsLvl[0], m_order);
    // USE updatePersistentCommunicationsLvl
  // if (Ncpu > 1) { m_mesh->communicationsPrimitives(m_eos, 0); }
  
  // //8) AMR initialization
  // //---------------------
  // m_mesh->procedureRaffinementInitialization(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_model, m_nbCellsTotalAMR, domains, m_eos, m_resumeSimulation, m_order);


std::cout<<"cpu "<<rankCpu //KS//BD//
<< " after balance"
<<std::endl;
for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
std::cout<<"cpu "<<rankCpu
<< " cellsLvl["<<lvl<<"].size() "<<cellsLvl[lvl].size()
<< " cellInterfacesLvl["<<lvl<<"].size() "<<cellInterfacesLvl[lvl].size()
<<std::endl;
}
  std::cout<<"cpu "<<rankCpu
    << " GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOD!!! "
    <<std::endl;
  MPI_Barrier(MPI_COMM_WORLD); //KS//BD//


  //Refinement and physical values to communicate

  //Delete cells in temporary cell buffer

  //NEED to verify final load...



  //END OF LOOP

  delete[] loadPerCPU;
  delete[] startPerCPU;
}

//***********************************************************************
