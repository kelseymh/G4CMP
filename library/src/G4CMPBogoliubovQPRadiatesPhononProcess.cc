/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRadiatesPhononProcess.cc
/// \brief Implementation of the G4CMPBogoliubovQPRadiatesPhononProcess class
//
// $Id$
//

#include "G4CMPBogoliubovQPRadiatesPhononProcess.hh"
#include "G4CMPBogoliubovQPRadiatesPhononRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4BogoliubovQP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor and destructor 
G4CMPBogoliubovQPRadiatesPhononProcess::G4CMPBogoliubovQPRadiatesPhononProcess(const G4String& aName)
  : G4VBogoliubovQPProcess(aName,fQPRadiatesPhononProcess)
{
  G4cout << "REL -- In QPRadiatesPhononProcess::Constructor()" << G4endl;
  UseRateModel(new G4CMPBogoliubovQPRadiatesPhononRate);
}

G4CMPBogoliubovQPRadiatesPhononProcess::~G4CMPBogoliubovQPRadiatesPhononProcess()
{
  G4cout << "REL -- In QPRadiatesPhononProcess::Destructor()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPBogoliubovQPRadiatesPhononProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPBogoliubovQPRadiatesPhononProcess::PostStepDoIt(const G4Track& aTrack,
								       const G4Step& aStep)
{
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRadiatesPhononProcess::PostStepDoIt ----------" << G4endl;
    G4cout << "PSDI Function Point A | poststeppoint velocity in QPRadiates poststepdoit is: " << aStep.GetPostStepPoint()->GetVelocity() << G4endl;
    G4cout << "PSDI Function Point A | track velocity in QPRadiates poststepdoit is: " << aTrack.GetVelocity() << G4endl;
  }
  

  aParticleChange.Initialize(aTrack);

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "PSDI Function Point B | poststeppoint velocity in QPRadiates poststepdoit, after initializing particle change, is: " << aStep.GetPostStepPoint()->GetVelocity() << G4endl;
    G4cout << "PSDI Function Point B | track velocity in QPRadiates poststepdoit, afteer initializing particle change, is: " << aTrack.GetVelocity() << G4endl;
  }

  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event -- if the code is working properly, this should
  //   basically never happen... REL May need to revisit this. I think it's still true but will have to check against the QP transport stuff.
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the Bogoliubov RadiatesPhonon process and we find ourselves on a boundary. Should this ever happen?";
    G4Exception("G4CMPBogoliubovQPRadiatesPhononProcess::PostStepDoIt", "BogoliubovQPRadiatesPhonon001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }

  
  
  //2. Identify the current QP's energy and velocity and use it to draw an energy for a radiated phonon.
  G4double qpEnergy = GetKineticEnergy(aTrack);
  G4double velocity = aStep.GetPostStepPoint()->GetVelocity(); 
  G4ThreeVector momDir = aTrack.GetMomentumDirection();
  
  //3. Actually generate that radiated phonon, and reduce the QP's energy.
  G4double radiatedPhonEnergy = PhononEnergyRand(qpEnergy);
  GenerateRadiatedPhonon(radiatedPhonEnergy,aTrack,aStep);
  aParticleChange.ProposeEnergy((qpEnergy-radiatedPhonEnergy));

  //4. Now we do the artificial setting of the QP's velocity and momentum direction again. These lines are aphysical but
  //   are okay because we are going to have to do the QP diffusion modeling in a hacky way anyway...
  aParticleChange.ProposeVelocity(velocity); 
  
  //  aParticleChange.ProposeMomentumDirection(momDir);
  RandomizeFinalStateMomentumDirectionInXY();
  
  //that is inconsistent with the diffusion physics. REL

  //4. Do the clear interaction lengths thing because we do still have a particle here.
  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  
  //5. Return the particle change
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4CMPBogoliubovQPRadiatesPhononProcess::IsApplicable(const G4ParticleDefinition& aPD) {

  // Allow all phonon types, because type is changed during tracking
  return G4VBogoliubovQPProcess::IsApplicable(aPD);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Take the phonon energy and produce a recombination phonon from this.
void G4CMPBogoliubovQPRadiatesPhononProcess::GenerateRadiatedPhonon(G4double phonEnergy,
								    const G4Track& aTrack,
								    const G4Step& aStep)
{
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRadiatesPhononProcess::GenerateRadiatedPhonon ----------" << G4endl;
  }
  
  //Now create the phonon
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(), theLattice->GetSTDOS(),theLattice->GetFTDOS());    
  G4ThreeVector dir1 = G4RandomDirection();    
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack,mode,dir1,phonEnergy,aTrack.GetGlobalTime(),aTrack.GetPosition());

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "GRP Function Point A | energy of radiated phonon: " << phonEnergy << ", whereas twice the gap is: " << 2.0*fGapEnergy << G4endl;
  }
  
  //Check to make sure the secondary was actually produced
  if (!sec1) {
    G4Exception("G4CMPBogoliubovQPRadiatesPhononProcess::GenerateRadiatedPhonon", "BogoliubovRadiatesPhonon002",
		JustWarning, "Error creating secondaries");
    return;
  }
  
  //Set the number of secondaries to 1 and add the secondary. Kill the QP.
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec1);
}
									
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Compute phonon energy distribution from quasiparticle in superconductor.
// This is the same as the original KaplanQP
G4double G4CMPBogoliubovQPRadiatesPhononProcess::PhononEnergyRand(G4double Energy) const
{
  
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-fGapEnergy*fGapEnergy/Energy))
  //           /
  //           sqrt((E'*E' - fGapEnergy*fGapEnergy);
  
  // Add buffer so first bin doesn't give zero denominator in pdfSum
  
  const G4double BUFF = 10000000000.; //REL used to be 1000 by default
  G4double xmin = fGapEnergy + fGapEnergy / BUFF;
  G4double xmax = Energy;
  G4double ymax = PhononEnergyPDF(Energy, xmin);
  
  G4double xtest = 0., ytest = ymax;
  do {
    ytest = G4UniformRand() * ymax;
    xtest = G4UniformRand() * (xmax - xmin) + xmin;
  } while (ytest > PhononEnergyPDF(Energy, xtest));  
  return Energy - xtest;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Phonon energy "pdf"
// This is the same as the original KaplanQP, but with a modification to fix a typo and
// bring this into agreement with the Kaplan paper
G4double G4CMPBogoliubovQPRadiatesPhononProcess::PhononEnergyPDF(G4double E, G4double x) const
{
    const G4double gapsq = fGapEnergy * fGapEnergy;
    return (x * (E - x) * (E - x) * (1 - gapsq / x / E) / sqrt(x * x - gapsq));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Pass-through to G4CMPVProcess class
G4double G4CMPBogoliubovQPRadiatesPhononProcess::GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond)
{
  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "---------- G4CMPBogoliubovQPRadiatesPhononProcess::GetMeanFreePath ----------" << G4endl;
  }
  
  //Need this to come first, so that it actually attempts a superconductor update.
  G4double mfpBase = G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);

  //Debugging
  if( verboseLevel > 5 ){
    G4cout << "GMFP Function Point A | Mean free path in QPRadiatesPhononProcess: " << mfpBase << ", with nMFPsLeft: " << GetNumberOfInteractionLengthLeft() << G4endl;
  }
  
  //If we don't trigger that exception, continue.
  return mfpBase;
}
