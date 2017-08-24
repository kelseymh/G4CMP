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

#include "G4PhononScattering.hh"
#include "G4CMPPhononScatteringRate.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"


G4PhononScattering::G4PhononScattering(const G4String& aName)
  : G4VPhononProcess(aName, fPhononScattering) {
  UseRateModel(new G4CMPPhononScatteringRate);
}

G4PhononScattering::~G4PhononScattering() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononScattering::PostStepDoIt( const G4Track& aTrack,
						     const G4Step& aStep) {
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
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

  //Initialize particle change
  aParticleChange.Initialize(aTrack);
  
  //randomly generate a new direction and polarization state
  G4ThreeVector newDir = G4RandomDirection();
  G4int polarization = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
					  theLattice->GetSTDOS(),
					  theLattice->GetFTDOS());

  if (verboseLevel>1) {
    G4cout << " Changing to "
	   << G4PhononPolarization::Get(polarization)->GetParticleName() << " "
	   << " toward " << newDir << G4endl;
  }

  // Generate the new track after scattering
  // FIXME:  If polarization state is the same, just step the track!
  if (verboseLevel) {
    G4cout << " Creating secondary using touchable for "
	   << aTrack.GetTouchable()->GetVolume()->GetName() << G4endl;
  }

  G4Track* sec =
    G4CMP::CreatePhonon(aTrack.GetTouchable(), polarization, newDir,
                        GetKineticEnergy(aTrack), aTrack.GetGlobalTime(),
                        aTrack.GetPosition());
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec);

  // Scattered phonon replaces current track
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
