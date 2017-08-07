/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPLukeEmissionRate.cc
/// \brief Compute emission rate for Luke-Neganov phonons.
//
// $Id$

#include "G4CMPLukeEmissionRate.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include <math.h>


// Scattering rate is computed from electric field

G4double G4CMPLukeEmissionRate::Rate(const G4Track& aTrack) const {
  const_cast<G4CMPLukeEmissionRate*>(this)->LoadDataForTrack(&aTrack);

  // Sanity check -- IsApplicable() should protect against this
  if (!G4CMP::IsChargeCarrier(aTrack)) {
    G4Exception("G4CMPLukeEmissionRate::Rate", "Luke001", EventMustBeAborted, 
		("Invalid particle "+aTrack.GetDefinition()->GetParticleName()).c_str());
    return 0.;
  }

  G4double kmag = 0.; G4double l0 = 0.; G4double mass = 0.;
  if (G4CMP::IsElectron(aTrack)) {
    kmag = theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
				     GetLocalVelocityVector(aTrack)).mag();
    l0 = theLattice->GetElectronScatter();
    mass = theLattice->GetElectronMass();	// Scalar mass
  } else if (G4CMP::IsHole(aTrack)) {
    kmag = GetLocalWaveVector(aTrack).mag();
    l0 = theLattice->GetHoleScatter();
    mass = theLattice->GetHoleMass();
  }

  if (verboseLevel > 1) 
    G4cout << "LukeEmissionRate kmag = " << kmag*m << " /m" << G4endl;

  G4double kSound = theLattice->GetSoundSpeed() * mass / hbar_Planck;

  if (kmag <= kSound) return 0.;
  
  // Time step corresponding to Mach number (avg. time between radiations)
  return 1./ChargeCarrierTimeStep(kmag/kSound, l0);
}
