/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononScattering.cc
/// \brief Implementation of the G4PhononScattering class
//
// $Id$
//
// 20131111  Add verbose output for MFP calculation
// 20140312  Follow name change CreateSecondary -> CreatePhonon
// 20140331  Add required process subtype code
// 20170620  Follow interface changes in G4CMPSecondaryUtils
// 20170805  Move GetMeanFreePath() to scattering-rate model
// 20170819  Overwrite track's particle definition instead of killing

#include "G4PhononScattering.hh"
#include "G4CMPPhononScatteringRate.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "Randomize.hh"


G4PhononScattering::G4PhononScattering(const G4String& aName)
  : G4VPhononProcess(aName, fPhononScattering) {
  UseRateModel(new G4CMPPhononScatteringRate);
}

G4PhononScattering::~G4PhononScattering() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononScattering::PostStepDoIt( const G4Track& aTrack,
						     const G4Step& aStep) {
  // Initialize particle change
  aParticleChange.Initialize(aTrack);
  
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary ||
      postStepPoint->GetStepStatus()==fWorldBoundary) {
    return &aParticleChange;			// Don't want to reset IL
  }

  if (verboseLevel) G4cout << GetProcessName() << "::PostStepDoIt" << G4endl;
  if (verboseLevel>1) {
    G4StepPoint* preStepPoint = aStep.GetPreStepPoint();
    G4cout << " Track " << aTrack.GetDefinition()->GetParticleName()
	   << " vol " << aTrack.GetTouchable()->GetVolume()->GetName()
	   << " prePV " << preStepPoint->GetPhysicalVolume()->GetName()
	   << " postPV " << postStepPoint->GetPhysicalVolume()->GetName()
	   << " step-length " << aStep.GetStepLength()
	   << G4endl;
  }

  // Randomly generate a new direction and polarization state
  G4ThreeVector newK = G4RandomDirection();
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
					       theLattice->GetSTDOS(),
					       theLattice->GetFTDOS());

  if (verboseLevel>1) {
    G4cout << " Changing to "
	   << G4PhononPolarization::Get(mode)->GetParticleName() << " "
	   << " toward " << newK << G4endl;
  }

  // Replace track's particle type according to new polarization
  if (mode != G4PhononPolarization::Get(aTrack.GetParticleDefinition())) {
    const G4ParticleDefinition* newPD = G4PhononPolarization::Get(mode);
    auto theDP = const_cast<G4DynamicParticle*>(aTrack.GetDynamicParticle());
    theDP->SetDefinition(newPD);

    if (verboseLevel>1) {		// Sanity check, report back PD
      G4cout << " track now " << aTrack.GetDefinition()->GetParticleName()
	     << G4endl;
    }
  }

  // Assign new wave vector direction to track (ought to happen later!)
  auto trkInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack);
  trkInfo->SetWaveVector(newK);

  // Set velocity and direction according to new wave vector direction
  G4double vgrp = theLattice->MapKtoV(mode, newK);
  G4ThreeVector vdir = theLattice->MapKtoVDir(mode, newK);
  RotateToGlobalDirection(vdir);

  if (verboseLevel>1)
    G4cout << " new vgrp " << vgrp << " along " << vdir << G4endl;

  aParticleChange.ProposeMomentumDirection(vdir);
  aParticleChange.ProposeVelocity(vgrp);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
