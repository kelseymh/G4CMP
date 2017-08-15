/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPDownconversionRate.cc
/// \brief Compute rate for phonon anharmonic decay (downconversion)
//
// $Id$
//
// 20170815  Drop call to LoadDataForTrack(); now handled in process.

#include "G4CMPDownconversionRate.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"


// Scattering rate is computed from electric field

G4double G4CMPDownconversionRate::Rate(const G4Track& aTrack) const {
  G4double A = theLattice->GetAnhDecConstant();
  G4double Eoverh = GetKineticEnergy(aTrack)/h_Planck;
  
  return (Eoverh*Eoverh*Eoverh*Eoverh*Eoverh*A);
}
