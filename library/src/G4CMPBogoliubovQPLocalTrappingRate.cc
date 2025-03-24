/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPLocalTrappingRate.cc
/// \brief Compute rate for QP local trapping 
//

#include "G4CMPBogoliubovQPLocalTrappingRate.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"



G4double G4CMPBogoliubovQPLocalTrappingRate::Rate(const G4Track& aTrack) const {

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPLocalTrappingRate::Rate ----------" << G4endl;
  }
  
  G4double localTrapTau = theLattice->GetSCQPLocalTrappingTau();

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "R Function Point A | qp local trapping tau: " << localTrapTau << G4endl;
  }

  //Some safeguards
  G4double rate = 0;
  if( localTrapTau <= 0.0 ){ rate = 0; }
  else{ rate = 1.0 / localTrapTau; }

  //More debugging
  if( verboseLevel > 5 ){
    G4cout << "R Function Point B | qp local trapping ultimate rate: " << rate << G4endl;
  }
  return rate;
}
