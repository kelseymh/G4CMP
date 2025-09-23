/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPRadiatesPhononProcess.cc
/// \brief Implementation of the G4CMPQPRadiatesPhononProcess class
//
// $Id$
//

#include "G4CMPQPRadiatesPhononProcess.hh"
#include "G4CMPQPRadiatesPhononRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4QP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"

// Constructor and destructor 
G4CMPQPRadiatesPhononProcess::
G4CMPQPRadiatesPhononProcess(const G4String& aName)
  : G4VQPProcess(aName,fQPRadiatesPhononProcess) {
  UseRateModel(new G4CMPQPRadiatesPhononRate);
}

G4CMPQPRadiatesPhononProcess::~G4CMPQPRadiatesPhononProcess() {
}


void G4CMPQPRadiatesPhononProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


G4VParticleChange*
G4CMPQPRadiatesPhononProcess::PostStepDoIt(const G4Track& aTrack,
					   const G4Step& aStep) {
  
  //Debugging
  if (verboseLevel > 5){
    G4cout << "-- G4CMPQPRadiatesPhononProcess::PostStepDoIt --" << G4endl;
    G4cout << "PSDI Function Point A | poststeppoint velocity in QPRadiates "
	   << "poststepdoit is: " << aStep.GetPostStepPoint()->GetVelocity()
	   << G4endl;
    G4cout << "PSDI Function Point A | track velocity in QPRadiates "
	   << "poststepdoit is: " << aTrack.GetVelocity() << G4endl;
  }
  
  aParticleChange.Initialize(aTrack);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "PSDI Function Point B | poststeppoint velocity in QPRadiates "
	   << "poststepdoit, after initializing particle change, is: "
	   << aStep.GetPostStepPoint()->GetVelocity() << G4endl;
    G4cout << "PSDI Function Point B | track velocity in QPRadiates "
	   << "poststepdoit, afteer initializing particle change, is: "
	   << aTrack.GetVelocity() << G4endl;
  }
  
  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event --
  //   if the code is working properly, this should basically never happen...
  //   REL May need to revisit this. I think it's still true but will have to
  //   check against the QP transport stuff.
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the Bogoliubov "
	<< "RadiatesPhonon process and we find ourselves on a boundary. Should "
	<< "this ever happen?";
    G4Exception("G4CMPQPRadiatesPhononProcess::PostStepDoIt",
		"QPRadiatesPhonon001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }
  
   
  //2. Identify the current QP's energy and velocity and use it to draw an
  //   energy for a radiated phonon.
  G4double qpEnergy = GetKineticEnergy(aTrack);
  G4double velocity = aStep.GetPostStepPoint()->GetVelocity(); 
  G4ThreeVector momDir = aTrack.GetMomentumDirection();
  
  //3. Actually generate that radiated phonon, and reduce the QP's energy.
  G4double radiatedPhonEnergy = PhononEnergyRand(qpEnergy);
  GenerateRadiatedPhonon(radiatedPhonEnergy,aTrack,aStep);
  aParticleChange.ProposeEnergy((qpEnergy-radiatedPhonEnergy));

  //4. Now we do the artificial setting of the QP's velocity and momentum
  //   direction again. Velocity should stay the same as it was.
  aParticleChange.ProposeVelocity(velocity);   
  RandomizeFinalStateMomentumDirectionInXY();
  
  //4. Do the clear interaction lengths thing because we do still have a
  //   particle here.
  ClearNumberOfInteractionLengthLeft();	// All processes should do this!

  //5. Return the particle change
  return &aParticleChange;
}

G4bool
G4CMPQPRadiatesPhononProcess::IsApplicable(const G4ParticleDefinition& aPD) {

  // Allow all phonon types, because type is changed during tracking
  return G4VQPProcess::IsApplicable(aPD);
}

// Take the phonon energy and produce a recombination phonon from this.
void G4CMPQPRadiatesPhononProcess::GenerateRadiatedPhonon(G4double phonEnergy,
							  const G4Track& aTrack,
							  const G4Step& aStep) {

  //Debugging
  if(verboseLevel > 5) {
    G4cout << "-- G4CMPQPRadiatesPhononProcess::GenerateRadiatedPhonon --"
	   << G4endl;
  }
  
  //Now create the phonon
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
					       theLattice->GetSTDOS(),
					       theLattice->GetFTDOS());    
  G4ThreeVector dir1 = G4RandomDirection();    
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack,mode,dir1,phonEnergy,
				      aTrack.GetGlobalTime(),
				      aTrack.GetPosition());

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "GRP Function Point A | energy of radiated phonon: "
	   << phonEnergy << ", whereas twice the gap is: " << 2.0*fGapEnergy
	   << G4endl;
  }
  
  //Check to make sure the secondary was actually produced
  if (!sec1) {
    G4Exception("G4CMPQPRadiatesPhononProcess::GenerateRadiatedPhonon",
		"QPRadiatesPhonon002", JustWarning,
		"Error creating secondaries");
    return;
  }
  
  //Set the number of secondaries to 1 and add the secondary. Kill the QP.
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec1);
}
									
// Compute phonon energy distribution from quasiparticle in superconductor.
// This is the same as the original KaplanQP
G4double G4CMPQPRadiatesPhononProcess::PhononEnergyRand(G4double Energy) const {
  
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

// Phonon energy "pdf"
// This is the same as the original KaplanQP, but with a modification to fix a
// typo and bring this into agreement with the Kaplan paper
G4double
G4CMPQPRadiatesPhononProcess::PhononEnergyPDF(G4double E, G4double x) const {
    const G4double gapsq = fGapEnergy * fGapEnergy;
    return (x * (E - x) * (E - x) * (1 - gapsq / x / E) / sqrt(x * x - gapsq));
}


//Pass-through to G4CMPVProcess class
G4double G4CMPQPRadiatesPhononProcess::GetMeanFreePath(const G4Track& trk,
						       G4double prevstep,
						       G4ForceCondition* cond) {

  //Debugging
  if(verboseLevel > 5) {
    G4cout << "-- G4CMPQPRadiatesPhononProcess::GetMeanFreePath --" << G4endl;
  }
  
  //Need this to come first, so that it actually attempts a superconductor
  //update.
  G4double mfpBase = G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "GMFP Function Point A | Mean free path in "
	   << "QPRadiatesPhononProcess: " << mfpBase << ", with nMFPsLeft: "
	   << GetNumberOfInteractionLengthLeft() << G4endl;
  }


  //If we don't trigger that exception, continue.
  return mfpBase;
}
