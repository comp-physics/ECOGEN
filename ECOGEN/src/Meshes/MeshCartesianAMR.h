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

#ifndef MESHCARTESIANAMR_H
#define MESHCARTESIANAMR_H

//! \file      MeshCartesian.h
//! \author    K. Schmidmayer, F. Petitpas, B. Dorschner
//! \version   1.0
//! \date      February 19 2019

#include "MeshCartesian.h"

class MeshCartesianAMR : public MeshCartesian
{
public:
  MeshCartesianAMR(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
    std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ,
		int lvlMax = 0, double criteriaVar = 1.e10, bool varRho = false, bool varP = false, bool varU = false, 
    bool varAlpha = false, double xiSplit = 1., double xiJoin = 1.);
  virtual ~MeshCartesianAMR();

  virtual int initializeGeometrie(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, bool pretraitementParallele, std::string ordreCalcul);
  void initializeGeometrieAMR(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string ordreCalcul);
  void createCellInterfacesFacesAndGhostCells(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string ordreCalcul, decomposition::Decomposition* _decomp);
	virtual void genereTableauxCellsCellInterfacesLvl(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::vector<Cell *> **cellsLvl,
		std::vector<CellInterface *> **cellInterfacesLvl);
  virtual void procedureRaffinementInitialization(std::vector<Cell *> *cellsLvl, std::vector<CellInterface *> *cellInterfacesLvl,
		const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, std::vector<GeometricalDomain*> &domains,
		TypeMeshContainer<Cell *> &cells, Eos **eos, const int &resumeSimulation);
  virtual void procedureRaffinement(std::vector<Cell *> *cellsLvl, std::vector<CellInterface *> *cellInterfacesLvl, const int &lvl,
    const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, TypeMeshContainer<Cell *> &cells, Eos **eos);
  virtual std::string whoAmI() const;

  //Printing / Reading
  virtual void ecritHeaderPiece(std::ofstream &fileStream, std::vector<Cell *> *cellsLvl) const;
  virtual void recupereNoeuds(std::vector<double> &jeuDonnees) const;
  virtual void recupereConnectivite(std::vector<double> &jeuDonnees) const;
  virtual void recupereOffsets(std::vector<double> &jeuDonnees) const;
  virtual void recupereTypeCell(std::vector<double> &jeuDonnees) const;
  virtual void recupereDonnees(std::vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase) const;
  virtual void setDataSet(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl, const int var, int phase) const;
  virtual void refineCell(Cell *cell, const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR);

  //Accesseurs
  virtual int getLvlMax() const { return m_lvlMax; };

	//Pour parallele
	virtual void initializePersistentCommunications(const int numberPhases, const int numberTransports, const TypeMeshContainer<Cell *> &cells, std::string ordreCalcul);
	virtual void communicationsPrimitives(const TypeMeshContainer<Cell *> &cells, Eos **eos, const int &lvl, Prim type = vecPhases);
	virtual void communicationsSlopes(const TypeMeshContainer<Cell *> &cells, const int &lvl);
	virtual void communicationsVector(const TypeMeshContainer<Cell *> &cells, std::string nameVector, const int &dim, const int &lvl, int num, int index);
	virtual void communicationsAddPhys(const std::vector<AddPhys*> &addPhys, const TypeMeshContainer<Cell *> &cells, const int &lvl);
  virtual void communicationsTransports(const TypeMeshContainer<Cell *> &cells, const int &lvl);
	virtual void finalizeParallele(const int &lvlMax);

private:
  int m_lvlMax;                              //!<Niveau maximal sur l arbre AMR (si m_lvlMax = 0, pas d AMR)
	double m_criteriaVar;                      //!<Valeur du criteria a depasser sur la variation d'une variable pour le (de)raffinement (met xi=1.)
	bool m_varRho, m_varP, m_varU, m_varAlpha; //!<Choix sur quelle variation on (de)raffine
	double m_xiSplit, m_xiJoin;                //!<Valeur de xi pour split ou join les mailles
  std::vector<Cell *> **m_cellsLvl;          //!<Pointer vers le tableau de vecteurs contenant les cells de compute, un vecteur par niveau.
	std::vector<Cell *> *m_cellsLvlGhost;      //!<Tableau de vecteurs contenant les cells fantomes, un vecteur par niveau.

};

#endif // MESHCARTESIANAMR_H