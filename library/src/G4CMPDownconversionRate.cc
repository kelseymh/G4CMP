/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPDownconversionRate.cc
/// \brief Compute rate for phonon anharmonic decay (downconversion)
//
// $Id$

#include "G4CMPDownconversionRate.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>


// Scattering rate is computed from electric field

G4double G4CMPDownconversionRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPDownconversionRate*>(this)->LoadDataForTrack(&aTrack);

  G4double A = theLattice->GetAnhDecConstant();
  G4double Eoverh = GetKineticEnergy(aTrack)/h_Planck;
  
  return (Eoverh*Eoverh*Eoverh*Eoverh*Eoverh*A);
}
