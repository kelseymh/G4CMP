/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPPhononScatteringRate.cc
/// \brief Compute rate for phonon impurity scattering (mode mixing)
//
// $Id$

#include "G4CMPPhononPolycrystalElasticScatteringRate.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"


// Scattering rate is computed from electric field

G4double G4CMPPhononPolycrystalElasticScatteringRate::Rate(const G4Track& aTrack) const {
  G4double scatMeanFreePath = theLattice->GetPolycrystalElasticScatterMFP();
  G4double scatMeanFreePath2 = theLattice->GetScatteringConstant();

  G4cout << "Scat mean free path: " << scatMeanFreePath << ", while scattering constant: " << scatMeanFreePath2 << G4endl; //REL probably should remove reference to scattering constant here (among other things)
  G4cout << "LDOS: " << theLattice->GetLDOS() << G4endl;
  
  //Some safeguards for cases where we don't want this to trigger (because it'll
  //try to run even for single crystals.
  G4double vel = aTrack.GetVelocity();
  G4cout << "vel: " << vel << G4endl;
  G4double rate = 0;
  if( scatMeanFreePath <= 0.0 ){ rate = 0; }
  else{ rate = vel / scatMeanFreePath; }

  G4cout << "ultimate rate: " << rate << G4endl;
  return rate;
}
