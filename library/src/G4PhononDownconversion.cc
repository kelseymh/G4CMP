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

#include "G4PhononDownconversion.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPDownconversionRate.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4CMPAnharmonicDecay.hh"
#include <cmath>


G4PhononDownconversion::G4PhononDownconversion(const G4String& aName)
  : G4VPhononProcess(aName, fPhononDownconversion),
    fBeta(0.), fGamma(0.), fLambda(0.), fMu(0.), fvLvT(1.) {
  UseRateModel(new G4CMPDownconversionRate);

#ifdef G4CMP_DEBUG
  if (verboseLevel) {
    output.open("phonon_downsampling_stats", std::ios_base::app);
    if (output.good()) {
      output << "First Daughter Theta,Second Daughter Theta,First Daughter Energy [eV],Second Daughter Energy [eV],"
	"Decay Branch,First Daughter Weight,Second Daughter Weight,Parent Weight,"
	"Number of Outgoing Tracks,Parent Energy [eV]\n";
    } else {
      G4cerr << "Could not open phonon debugging output file!" << G4endl;
    }
  }
#endif
}

G4PhononDownconversion::~G4PhononDownconversion() {
#ifdef G4CMP_DEBUG
  output.close();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononDownconversion::PostStepDoIt( const G4Track& aTrack,
							 const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
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

  // Only longitudinal phonons decay
  if (aTrack.GetDefinition() != G4PhononLong::Definition()) {
    return &aParticleChange;		// Don't reset interaction length!
  }

  G4CMPAnharmonicDecay::DoDecay(aTrack, aStep, aParticleChange);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PhononDownconversion::IsApplicable(const G4ParticleDefinition& aPD) {
  // Only L-phonons decay
  /***** , but need to check actively changing phonon type
  return (&aPD==G4PhononLong::PhononDefinition());
  *****/
  return G4VPhononProcess::IsApplicable(aPD);
}
