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

#include "G4PhononScattering.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononLong.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"


G4PhononScattering::G4PhononScattering(const G4String& aName)
  : G4VPhononProcess(aName, fPhononScattering) {;}

G4PhononScattering::~G4PhononScattering() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PhononScattering::GetMeanFreePath(const G4Track& aTrack,
					     G4double /*previousStepSize*/,
					     G4ForceCondition* condition) {
  //Dynamical constants retrieved from PhysicalLattice
  G4double B = theLattice->GetScatteringConstant();
  G4double Eoverh = aTrack.GetKineticEnergy()/h_Planck;

  //Calculate mean free path
  G4double mfp = aTrack.GetVelocity()/(Eoverh*Eoverh*Eoverh*Eoverh*B);

  if (verboseLevel > 1)
    G4cout << "G4PhononScattering::GetMeanFreePath = " << mfp << G4endl;

  *condition = NotForced;
 
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononScattering::PostStepDoIt( const G4Track& aTrack,
						     const G4Step& aStep) {
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  
  //Initialize particle change
  aParticleChange.Initialize(aTrack);
  
  //randomly generate a new direction and polarization state
  G4ThreeVector newDir = G4RandomDirection();
  G4int polarization = ChoosePolarization(theLattice->GetLDOS(),
					  theLattice->GetSTDOS(),
					  theLattice->GetFTDOS());

  // Generate the new track after scattering
  // FIXME:  If polarization state is the same, just step the track!
  G4Track* sec =
    CreatePhonon(polarization, newDir, aTrack.GetKineticEnergy());
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec);

  // Scattered phonon replaces current track
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
