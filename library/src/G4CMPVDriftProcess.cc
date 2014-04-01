//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/src/G4CMPVDriftProcess.cc
/// \brief Implementation of the G4CMPVDriftProcess base class
//
// $Id$
//
// 20140325  Move time-step calculation here from TimeStepper and LukeScat
// 20140331  Add required process subtype code

#include "G4CMPVDriftProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessType.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"


// Constructor and destructor
// NOTE:  Initial values are arbitrary and non-physical
G4CMPVDriftProcess::G4CMPVDriftProcess(const G4String& processName,
				       G4CMPProcessSubType stype)
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils(),
    velLong(330*m/s), mc_e(electron_mass_c2), l0_e(1*nm), 
    ksound_e(velLong*mc_e/hbar_Planck), mc_h(mc_e), l0_h(l0_e),
    ksound_h(ksound_e) {
  SetProcessSubType(stype);
#ifdef G4CMP_DEBUG
  G4cout << GetProcessName() << " is created " << G4endl;
#endif
}

G4CMPVDriftProcess::~G4CMPVDriftProcess() {;}


// Only applies to the known charge carriers

G4bool G4CMPVDriftProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4CMPDriftElectron::Definition() ||
	  &aPD==G4CMPDriftHole::Definition() );
}


// Get additional parameters from lattice for carriers

void G4CMPVDriftProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);

  velLong = theLattice->GetSoundSpeed();

  mc_e = theLattice->GetElectronMass();		// Scalar mass (Henning-Vogt)
  l0_e = theLattice->GetElectronScatter();
  ksound_e = velLong * mc_e/hbar_Planck;	// Wavevector for e @ "Mach 1"

  mc_h = theLattice->GetHoleMass();
  l0_h = theLattice->GetHoleScatter();
  ksound_h = velLong * mc_h/hbar_Planck;	// Wavevector for h @ "Mach 1"
}


// Initialize wave vectors for currently active track(s)

void G4CMPVDriftProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  LoadDataForTrack(track);
}

void G4CMPVDriftProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  ReleaseTrack();
}


// Compute characteristic time step for charge carrier
// Parameters are "Mach number" (ratio with sound speed) and scattering length

G4double 
G4CMPVDriftProcess::ChargeCarrierTimeStep(G4double mach, G4double l0) const {
  if (mach < 1.) return 3*l0/velLong;	// Sanity check if below sound speed

  return ( 3*l0 * mach / (velLong * (mach-1)*(mach-1)*(mach-1)) );
}
