/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononDownconversion.cc
/// \brief Implementation of the G4PhononDownconversion class
//
// $Id$
//
// 20131111  Add verbose output for MFP calculation
// 20131115  Initialize data buffers in ctor
// 20140312  Follow name change CreateSecondary -> CreatePhonon
// 20140331  Add required process subtype code
// 20160624  Use GetTrackInfo() accessor
// 20161114  Use new PhononTrackInfo
// 20170620  Follow interface changes in G4CMPSecondaryUtils
// 20170801  Protect PostStepDoIt() from being called at boundary
// 20170802  Use G4CMP_DOWN_SAMPLE biasing with ChooseWeight(), move outside
//		of sub-functions.
// 20170805  Replace GetMeanFreePath() with scattering-rate model
// 20170820  Compute MFP for all phonon types, check for L-type in PostStep
// 20170821  Move hard-coded constants to lattice configuration
// 20170824  Add diagnostic output
// 20170928  Hide "output" usage behind verbosity check, as well as G4CMP_DEBUG
// 20181010  J. Singh -- Move functionality to G4CMPAnharmonicDecay.
// 20181011  M. Kelsey -- Add LoadDataForTrack() to initialize decay utility.
// 20191014  G4CMP-179:  Drop sampling of anharmonic decay (downconversion)
// 20200604  G4CMP-208:  Report accept-reject values of u,x,q for debugging.
// 20201109  Move debugging output creation to PostStepDoIt to allows settting
//		process verbosity via macro commands.
// 20220712  M. Kelsey -- Pass process pointer to G4CMPAnharmonicDecay

#include "G4PhononDownconversion.hh"
#include "G4CMPAnharmonicDecay.hh"
#include "G4CMPDownconversionRate.hh"
#include "G4PhononLong.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Constructor and destructor 

G4PhononDownconversion::G4PhononDownconversion(const G4String& aName)
  : G4VPhononProcess(aName, fPhononDownconversion),
    anharmonicDecay(new G4CMPAnharmonicDecay(this)) {
  UseRateModel(new G4CMPDownconversionRate);
}

G4PhononDownconversion::~G4PhononDownconversion() {
  delete anharmonicDecay;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Configure for current track including AnharmonicDecay utility

void G4PhononDownconversion::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);
  anharmonicDecay->LoadDataForTrack(track);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Pass verbosity through to decay utility

void G4PhononDownconversion::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
  anharmonicDecay->SetVerboseLevel(vb);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononDownconversion::PostStepDoIt(const G4Track& aTrack,
							const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    return &aParticleChange;			// Don't want to reset IL
  }

  // Only longitudinal phonons decay
  if (aTrack.GetDefinition() != G4PhononLong::Definition()) {
    return &aParticleChange;		// Don't reset interaction length!
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

  anharmonicDecay->DoDecay(aTrack, aStep, aParticleChange);
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PhononDownconversion::IsApplicable(const G4ParticleDefinition& aPD) {
  // Allow all phonon types, because type is changed during tracking
  return G4VPhononProcess::IsApplicable(aPD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

