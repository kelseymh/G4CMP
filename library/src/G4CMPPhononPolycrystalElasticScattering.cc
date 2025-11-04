/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPPhononPolycrystalElasticScattering.cc
/// \brief Implementation of the G4CMPPhononPolycrystalElasticScattering class
//

#include "G4CMPPhononPolycrystalElasticScattering.hh"
#include "G4CMPPhononPolycrystalElasticScatteringRate.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticeManager.hh"
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


G4CMPPhononPolycrystalElasticScattering::G4CMPPhononPolycrystalElasticScattering(const G4String& aName)
  : G4VPhononProcess(aName, fPhononPolycrystalElasticScattering) {
  UseRateModel(new G4CMPPhononPolycrystalElasticScatteringRate);
}

G4CMPPhononPolycrystalElasticScattering::~G4CMPPhononPolycrystalElasticScattering() {
}


// Can basically keep the same scattering physics that is in the
// G4CMPPhononPolycrystalElasticScattering for now at least.
// May want to slightly modify this after having a discussion about how the
// polycrystalline shape affects the orientation of the caustics, mode mixing,
// etc. 
G4VParticleChange* G4CMPPhononPolycrystalElasticScattering::
PostStepDoIt( const G4Track& aTrack, const G4Step& aStep) {
  
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

//Pass-through to G4CMPVProcess class
G4double G4CMPPhononPolycrystalElasticScattering::
GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond) {
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}
