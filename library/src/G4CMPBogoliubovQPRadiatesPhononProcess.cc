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
#include "G4LatticePhysical.hh"
#include "G4CMPSecondaryUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor and destructor 
G4CMPBogoliubovQPRadiatesPhononProcess::G4CMPBogoliubovQPRadiatesPhononProcess(const G4String& aName)
  : G4VBogoliubovQPProcess(aName,fQPRadiatesPhononProcess)
{  
  UseRateModel(new G4CMPBogoliubovQPRadiatesPhononRate);
}

G4CMPBogoliubovQPRadiatesPhononProcess::~G4CMPBogoliubovQPRadiatesPhononProcess()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPBogoliubovQPRadiatesPhononProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPBogoliubovQPRadiatesPhononProcess::PostStepDoIt(const G4Track& aTrack,
								       const G4Step& aStep)
{
  G4cout << "REL in BogoliubovQPRadiatesPhonon process poststepdoit." << G4endl;
  
  aParticleChange.Initialize(aTrack);
  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, kill the event -- if the code is working properly, this should
  //   basically never happen...
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    G4ExceptionDescription msg;
    msg << "For some reason we're running post-step do it for the Bogoliubov RadiatesPhonon process and we find ourselves on a boundary. Since the QPs are artificially told to move attometers by design, I'm not expecting this to happen.";
    G4Exception("G4CMPBogoliubovQPRadiatesPhononProcess::PostStepDoIt", "BogoliubovQPRadiatesPhonon001",EventMustBeAborted,msg);
    return &aParticleChange;		
  }

  //2. Identify the current QP's energy and velocity and use it to draw an energy for a radiated phonon.
  G4double qpEnergy = GetKineticEnergy(aTrack);
  G4double velocity = aTrack.GetVelocity(); // note: need to convert to m/s
  G4ThreeVector momDir = aTrack.GetMomentumDirection();
  
  //3. Actually generate that radiated phonon, and reduce the QP's energy.
  G4double radiatedPhonEnergy = PhononEnergyRand(qpEnergy);
  GenerateRadiatedPhonon(radiatedPhonEnergy,aTrack,aStep);
  aParticleChange.ProposeEnergy((qpEnergy-radiatedPhonEnergy));

  //4. Now we do the artificial setting of the phonon's velocity and momentum direction again. These lines are aphysical but
  //   are okay because we are going to have to do the QP diffusion modeling in a hacky way anyway...
  aParticleChange.ProposeMomentumDirection(momDir);
  aParticleChange.ProposeVelocity(velocity); //Hopefully this still continues to not give problems...

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

  //Now create the phonon
  G4int mode = G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(), theLattice->GetSTDOS(),theLattice->GetFTDOS());    
  G4ThreeVector dir1 = G4RandomDirection();    
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack,mode,dir1,phonEnergy,aTrack.GetGlobalTime(),aTrack.GetPosition());

  G4cout << "REL energy of radiated phonon: " << phonEnergy << ", whereas twice the gap is: " << 2.0*fGapEnergy << G4endl;
  
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
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}
