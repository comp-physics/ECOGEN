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

//! \file      FluxMultiP.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      June 5 2017

#include <cmath>
#include <algorithm>
#include "FluxMultiP.h"
#include "../Mixture.h"

FluxMultiP *fluxBufferMultiP;
FluxMultiP *sourceConsMultiP;

//***********************************************************************

FluxMultiP::FluxMultiP(){}

//***********************************************************************

FluxMultiP::FluxMultiP(ModMultiP *model, const int &numberPhases) : m_model(model)
{
  m_alpha = new double[numberPhases];
  m_masse = new double[numberPhases];
  m_energ = new double[numberPhases];
}

//***********************************************************************

FluxMultiP::~FluxMultiP()
{
  delete[] m_alpha;
  delete[] m_masse;
  delete[] m_energ;
}

//***********************************************************************

void FluxMultiP::printFlux() const
{
  std::cout << m_masse << " " << m_qdm.getX() << " " << m_energ << std::endl;
}

//***********************************************************************

void FluxMultiP::addFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] += coefA*fluxBufferMultiP->m_alpha[k];
    m_masse[k] += coefA*fluxBufferMultiP->m_masse[k];
    m_energ[k] += coefA*fluxBufferMultiP->m_energ[k];
  }
  m_qdm += coefA*fluxBufferMultiP->m_qdm;
  m_energMixture += coefA*fluxBufferMultiP->m_energMixture;
}

//***********************************************************************

void FluxMultiP::subtractFlux(double coefA, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] -= coefA*fluxBufferMultiP->m_alpha[k];
    m_masse[k] -= coefA*fluxBufferMultiP->m_masse[k];
    m_energ[k] -= coefA*fluxBufferMultiP->m_energ[k];
  }
  m_qdm -= coefA*fluxBufferMultiP->m_qdm;
  m_energMixture -= coefA*fluxBufferMultiP->m_energMixture;
}

//***********************************************************************

void FluxMultiP::multiply(double scalar, const int &numberPhases)
{
    for(int k=0;k<numberPhases;k++)
    {
      m_alpha[k] *= scalar;
      m_masse[k] *= scalar;
      m_energ[k] *= scalar;
    }
    m_qdm  *= scalar;
    m_energMixture *= scalar;
}

//***********************************************************************

void FluxMultiP::setBufferFlux(Cell &cell, const int &numberPhases)
{
  fluxBufferMultiP->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxMultiP::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
	double energieInterne(0.);
  for (int k = 0; k < numberPhases; k++)
	{
		m_alpha[k] = phases[k]->getAlpha();
		m_masse[k] = phases[k]->getAlpha() * phases[k]->getDensity();
		energieInterne = phases[k]->getEos()->computeEnergy(phases[k]->getDensity(), phases[k]->getPressure());
		m_energ[k] = phases[k]->getAlpha() * phases[k]->getDensity() * energieInterne;
	}
	m_qdm = mixture->getDensity()*mixture->getVelocity();
  m_energMixture = mixture->getDensity()*mixture->getTotalEnergy();
}

//***********************************************************************

void FluxMultiP::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  double pressure(0.), energieInterne(0.), rhoMel(0.);

  //Verification and correction if needed (alpha and mass for order 2)
  double un(0.);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] < 0.) m_alpha[k] = 0.;
      if (m_alpha[k] > 1.) m_alpha[k] = 1.;
      un += m_alpha[k];
    }
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    for (int k = 0; k < numberPhases; k++) {
      if (m_alpha[k] <= 1.e-15) m_alpha[k] = 1e-15;
      if (m_alpha[k] >= 1.-1.e-15) m_alpha[k] = 1.0 - 1e-15;
      un += m_alpha[k];
    }
  }
  for (int k = 0; k < numberPhases; k++) { m_alpha[k] /= un; }

  //Phases and mixture variables
  for (int k = 0; k < numberPhases; k++) {
      rhoMel = rhoMel + m_masse[k];
      phases[k]->setAlpha(m_alpha[k]);
      phases[k]->setDensity(m_masse[k] / std::max(m_alpha[k], epsilonAlphaNull));
      //Calcul Pressure
      energieInterne = m_energ[k] / std::max(m_masse[k], epsilonAlphaNull);
      pressure = TB->eos[k]->computePressure(phases[k]->getDensity(), energieInterne);
      phases[k]->setPressure(pressure);
      phases[k]->verifyAndCorrectPhase();
  }
  mixture->setVelocity(m_qdm.getX() / rhoMel, m_qdm.getY() / rhoMel, m_qdm.getZ() / rhoMel);
  //Erasing small velocity variations
  if (std::fabs(mixture->getU()) < 1.e-8) mixture->setU(0.);
  if (std::fabs(mixture->getV()) < 1.e-8) mixture->setV(0.);
  if (std::fabs(mixture->getW()) < 1.e-8) mixture->setW(0.);
  for (int k = 0; k < numberPhases; k++) {
    phases[k]->extendedCalculusPhase(mixture->getVelocity());
  }
  mixture->computeMixtureVariables(phases, numberPhases);
  //Reconstruction of total energy from total energy equation
  double totalEnergyMixture(0.);
  totalEnergyMixture = m_energMixture / rhoMel;
  mixture->setTotalEnergy(totalEnergyMixture);
}

//***********************************************************************

void FluxMultiP::setToZero(const int &numberPhases)
{
  for(int k=0;k<numberPhases;k++){
    m_alpha[k] = 0.;
    m_masse[k] = 0.;
    m_energ[k] = 0.;
  }
  m_qdm = 0.;
  m_energMixture = 0.;
}

//***********************************************************************

void FluxMultiP::setToZeroBufferFlux(const int &numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferMultiP->m_alpha[k] = 0.;
    fluxBufferMultiP->m_masse[k] = 0.;
    fluxBufferMultiP->m_energ[k] = 0.;
  }
  fluxBufferMultiP->m_qdm = 0.;
  fluxBufferMultiP->m_energMixture = 0.;
}

//***********************************************************************

void FluxMultiP::addNonCons(double coefA, const Cell *cell, const int &numberPhases)
{
  Phase *phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] += -coefA*phase->getAlpha()*fluxBufferMultiP->m_sM;
    m_energ[k] += coefA*phase->getAlpha()*phase->getPressure()*fluxBufferMultiP->m_sM;
  }
}

//***********************************************************************

void FluxMultiP::subtractNonCons(double coefA, const Cell *cell, const int &numberPhases)
{
  Phase *phase;
  for(int k=0;k<numberPhases;k++){
    phase = cell->getPhase(k);
    m_alpha[k] -= -coefA*phase->getAlpha()*fluxBufferMultiP->m_sM;
    m_energ[k] -= coefA*phase->getAlpha()*phase->getPressure()*fluxBufferMultiP->m_sM;
  }
}

//***********************************************************************

void FluxMultiP::schemeCorrection(Cell *cell, const int &numberPhases, Prim type) const
{
  double rhoN(0.), rhoNp1(0.), sommeDeltaek(0.);
  for (int k = 0; k < numberPhases; k++) {
    rhoN += fluxBufferMultiP->getMasse(k);
    rhoNp1 += m_masse[k];
    TB->Deltaek[k] = (m_energ[k] - fluxBufferMultiP->getEnergy(k));
    sommeDeltaek += TB->Deltaek[k];
  }

  double DeltaE, DeltarhoU2;
  DeltaE = m_energMixture - fluxBufferMultiP->getEnergyMix();
  DeltarhoU2 = m_qdm.squaredNorm() / rhoNp1 - fluxBufferMultiP->getQdm().squaredNorm() / rhoN;

  if (std::fabs(sommeDeltaek) > 1.e-10*m_energMixture) {
    for (int k = 0; k < numberPhases; k++) {
      m_energ[k] = fluxBufferMultiP->getEnergy(k) + TB->Deltaek[k] / sommeDeltaek * (DeltaE - 0.5*DeltarhoU2);
    }
  }
}

//***********************************************************************

void FluxMultiP::addSymmetricTerms(Phase **phases, Mixture *mixture, const int &numberPhases, const double &r, const double &v, const double &dt)
{
  for (int k = 0; k<numberPhases; k++)
  {
    //Note: there is no cylindrical or spherical terms in a transport equation
    m_masse[k] += -phases[k]->getAlpha() * phases[k]->getDensity() * v / r;
    m_energ[k] += -(phases[k]->getAlpha() * phases[k]->getDensity() * phases[k]->getEnergy() + phases[k]->getAlpha() * phases[k]->getPressure()) * v / r;
  }
  m_qdm += -v / r * mixture->getDensity() * mixture->getVelocity();
  m_energMixture += -(mixture->getDensity() * mixture->getTotalEnergy() + mixture->getPressure()) * v / r;
}

//***********************************************************************

void FluxMultiP::integrateSourceTermsGravity(Cell *cell, const double &dt, const int &numberPhases, const int &axe, const int &direction, const Coord &g)
{
  sourceConsMultiP->setToZero(numberPhases);
  //Mass and velocity extraction
  double rho = cell->getMixture()->getDensity();
  Coord u = cell->getMixture()->getVelocity();

  //Gravity force and work
  sourceConsMultiP->m_qdm = rho*g;
  sourceConsMultiP->m_energMixture = rho*Coord::scalarProduct(g, u);

  //Euler integration (order 1)
  m_qdm += dt*sourceConsMultiP->m_qdm;
  m_energMixture += dt*sourceConsMultiP->m_energMixture;
}

//***********************************************************************

void FluxMultiP::integrateSourceTermsHeating(Cell *cell, const double &dt, const int &numberPhases, const double &q)
{
  sourceConsMultiP->setToZero(numberPhases);

  //Version 1 on 2p model, then relax
  sourceConsMultiP->m_energMixture = q;
  double sumMasses(0.);
  for (int k = 0; k < numberPhases; k++) {
    sumMasses += m_masse[k];
  }
    for (int k = 0; k < numberPhases; k++) {
    sourceConsMultiP->m_energ[k] = m_masse[k] / sumMasses * q;
    m_energ[k] += dt*sourceConsMultiP->m_energ[k];
  }
  m_energMixture += dt*sourceConsMultiP->m_energMixture;


  ////Version 2 directly on Kapila model
  //double p0 = cell->getMixture()->getPressure();
  //double v = 1./cell->getMixture()->getDensity();
  //vector<double> vk(numberPhases);
  //vector<double> vk0(numberPhases);
  //for (int k = 0; k < numberPhases; k++) {
  //  vk0[k] = 1./cell->getPhase(k)->getDensity();
  //}
  //double pStar(p0+1.);
  //double dvkdp(0.);

  //int iteration(0);
  //double f(0.), df(1.);
  //do {
  //  pStar -= f / df; iteration++;
  //  if (iteration > 50) {
  //    errors.push_back(Errors("not converged in integrateSourceTermsHeating", __FILE__, __LINE__));
  //    break;
  //  }
  //  //Physical pressure?
  //  for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pStar); }
  //  //Specific volumes after heating
  //  f = 1.;
  //  df = 0.;
  //  for (int k = 0; k < numberPhases; k++) {
  //    vk[k] = TB->eos[k]->computeSpecificVolumeQ(p0, vk0[k], pStar, q*v*dt, &dvkdp);
  //    f -= m_masse[k] * vk[k];
  //    df -= m_masse[k] * dvkdp;
  //  }
  //} while (std::fabs(f)>1e-12);

  ////Cell update
  //Phase *phase;
  //double ak;
  //for (int k = 0; k < numberPhases; k++) {
  //  phase = cell->getPhase(k);
  //  ak = m_masse[k] * vk[k];
  //  phase->setAlpha(ak);
  //  phase->setDensity(1./vk[k]);
  //  phase->setPressure(pStar);

  //  //m_alpha[k] = m_masse[k] * vk[k];
  //  //m_energ[k] = m_masse[k] * TB->eos[k]->computeEnergy(1. / vk[k], pStar);
  //}

  //
  ////m_energMixture += dt*q;

  ////cell->getMixture()->setPressure(pStar);
  //cell->fulfillState();
}

//***********************************************************************

void FluxMultiP::integrateSourceTermsMRF(Cell *cell, const double &dt, const int &numberPhases, const Coord &omega)
{
  sourceConsMultiP->setToZero(numberPhases);
  //Mass and velocity extraction
  double rho = cell->getMixture()->getDensity();
  Coord u = cell->getMixture()->getVelocity();
  
  //Coriolis acceleration
  sourceConsMultiP->m_qdm = -2.*rho*Coord::crossProduct(omega,u);
  //Centrifugal acceleration
  sourceConsMultiP->m_qdm -= rho*Coord::crossProduct(omega, Coord::crossProduct(omega, cell->getPosition()));
  //Centrifugal acceleration work
  sourceConsMultiP->m_energMixture = Coord::scalarProduct(u, sourceConsMultiP->m_qdm);

  //Euler integration (order 1)
  m_qdm += dt*sourceConsMultiP->m_qdm;
  m_energMixture += dt*sourceConsMultiP->m_energMixture;
}

//***********************************************************************

void FluxMultiP::setCons(const Flux *cons, const int &numberPhases)
{
  for (int k = 0; k < numberPhases; k++) {
    m_alpha[k] = cons->getAlpha(k);
    m_masse[k] = cons->getMasse(k);
    m_energ[k] = cons->getEnergy(k);
  }
  m_qdm = cons->getQdm();
  m_energMixture = cons->getEnergyMix();
}

//***********************************************************************