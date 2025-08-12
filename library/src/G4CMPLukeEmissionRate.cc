/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPLukeEmissionRate.cc
/// \brief Compute emission rate for Luke-Neganov phonons.
//
// $Id$
//
// 20170815  Drop call to LoadDataForTrack(); now handled in process.
// 20170913  Check for electric field; compute "rate" to get up to Vsound
// 20170917  Add interface for threshold identification
// 20240207  Replacing wave vector with speed to Emission calculations

#include "G4CMPLukeEmissionRate.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include <math.h>


// Scattering rate is computed from electric field

G4double G4CMPLukeEmissionRate::Rate(const G4Track& aTrack) const {
  // Sanity check -- IsApplicable() should protect against this
  if (!G4CMP::IsChargeCarrier(aTrack)) {
    G4Exception("G4CMPLukeEmissionRate::Rate", "Luke001", EventMustBeAborted, 
		("Invalid particle "+aTrack.GetDefinition()->GetParticleName()).c_str());
    return 0.;
  }

  G4double kSound = 0.; G4double mass = 0.; G4double l0 = 0.; G4double vsound = 0.;
  G4ThreeVector ktrk(0.);
  G4ThreeVector ptrk = GetLocalMomentum(aTrack);
  G4int iValley = GetValleyIndex(aTrack);
  if (G4CMP::IsElectron(aTrack)) {
    l0 = theLattice->GetElectronScatter();
    ktrk = theLattice->MapPtoK(iValley, ptrk);
    // Turning wavevector to spherical frame where electrons act like holes
    // as the mass is isotropic
    ktrk = theLattice->EllipsoidalToSphericalTranformation(iValley, ktrk);
    mass = sqrt(theLattice->GetElectronMass()*theLattice->GetElectronDOSMass());
    vsound = theLattice->GetAverageSoundSpeed();
    // The l0 in configuration file is calculated using the density of states mass
    // l0 = l0*pow(theLattice->GetElectronMass(),3)/(pow(mass,3));
  } else if (G4CMP::IsHole(aTrack)) {
    l0 = theLattice->GetHoleScatter();
    ktrk = GetLocalWaveVector(aTrack);
    mass = theLattice->GetHoleMass();
    vsound = theLattice->GetSoundSpeed();
  }
  G4double kmag = ktrk.mag();

  G4double gammaSound = 1/sqrt(1.-vsound*vsound/c_squared);
  kSound = gammaSound*vsound*mass/hbar_Planck;

  return (kmag > kSound) ? 1./ChargeCarrierTimeStep(kmag/kSound, l0, vsound) : 0.;
}


// Energy threshold occurs at sound speed in material

G4double G4CMPLukeEmissionRate::Threshold(G4double Eabove) const {
  const G4Track* trk = GetCurrentTrack();	// For convenience below

  G4ThreeVector vtrk = GetLocalVelocityVector(trk);
  G4double vsound = theLattice->GetSoundSpeed();
  G4double Esound = 0.; G4double mass = 0.;
  if (G4CMP::IsElectron(trk)) {
    mass = theLattice->GetElectronMass();
  } else {
    mass = theLattice->GetHoleMass();
  }

  // The sound speed is in the spherical frame, where mass is isotropic
  G4double gammaSound = 1/sqrt(1.-vsound*vsound/c_squared);
  Esound = (gammaSound-1)*mass*c_squared;

  if (verboseLevel>1) {
    G4cout << "G4CMPLukeEmissionRate::Threshold vtrk " << vtrk.mag()/(m/s)
	   << " vsound " << vsound/(m/s) << " m/s Esound " << Esound/eV
	   << " eV" << G4endl;
  }

  // Thresholds or pseudothresholds at multiples of Esound
  const G4double eStep = 25.;
  G4double ratio = Eabove/Esound;
  if (ratio > 1.) {
    ratio += (std::fmod(ratio,eStep)<0.95) ? 0. : eStep;  // Avoid Zeno paradox
    ratio = std::ceil(ratio/eStep) * eStep;

    if (verboseLevel>2) G4cout << " scaling by " << ratio << G4endl;
    
    Esound *= ratio;
  }

  return (Eabove < Esound) ? Esound : 0.;	// No thresholds above sound
}
