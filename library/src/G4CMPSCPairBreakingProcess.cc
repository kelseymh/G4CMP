/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPSCPairBreakingProcess.cc
/// \brief Implementation of the G4CMPSCPairBreakingClass
//
// $Id$
//

#include "G4CMPSCPairBreakingProcess.hh"
#include "G4CMPSCPairBreakingRate.hh"
#include "G4CMPSCUtils.hh"
#include "G4BogoliubovQP.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4RandomDirection.hh"
#include "G4CMPSecondaryUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor and destructor 
G4CMPSCPairBreakingProcess::G4CMPSCPairBreakingProcess(const G4String& aName)
  : G4VPhononProcess(aName,fSCPairBreakingProcess)
{  
  UseRateModel(new G4CMPSCPairBreakingRate);
  G4cout << "REL HereA_SCPairBreakingProcess" << G4endl;
}

G4CMPSCPairBreakingProcess::~G4CMPSCPairBreakingProcess()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPSCPairBreakingProcess::SetVerboseLevel(G4int vb) {
  verboseLevel = vb;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VParticleChange* G4CMPSCPairBreakingProcess::PostStepDoIt(const G4Track& aTrack,
							    const G4Step& aStep) {


  aParticleChange.Initialize(aTrack);
  
  //Pseudocode
  //1. Determine if we're on a boundary surface. If we are, just return the particlechange that's initialized already and don't reset
  //   the interaction lengths (interaction lengths are subtracted deep in the G4 code)
  G4StepPoint * postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus() == fGeomBoundary ||
      postStepPoint->GetStepStatus() == fWorldBoundary) {
    return &aParticleChange;		
  }
  
   
  //2. Using the phonon energy, compute the energy of the quasiparticles that are produced from our G4CMPSCUtils class, given
  //   the effective temperature we're at.
  double phononEnergy = aTrack.GetKineticEnergy();
  std::pair<G4double,G4double> QPenergies = FetchQPEnergies(phononEnergy);
  
  //3. Using the two above-computed energies, generate the two secondaries (G4BogoliubovQPs) we want.
  GenerateBogoliubovQPPair(QPenergies,aTrack,aStep);

  //4. Print debugging info
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

  //5. Sanity check to make sure we actually generated two secondaries, and then kill the track
  if (aParticleChange.GetNumberOfSecondaries() == 2) {
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }
  else{
    G4cout << "----> REL Uh oh. Bogoliubov secondaries not produced somehow?" << G4endl;
  }

  //6. Return the particle change
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4CMPSCPairBreakingProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  // Allow all phonon types, because type is changed during tracking
  return G4VPhononProcess::IsApplicable(aPD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//For a given phonon energy, sample QP energies from the relevant distribution. The
//distributions are stored in a G4SCUtils class so we only have to integrate once
//and then sample as we go.
std::pair<G4double,G4double> G4CMPSCPairBreakingProcess::FetchQPEnergies(G4double phonEnergy)
{
  //From the G4CMPSUtils class, sample a qp energy:
  G4double qp1Energy = QPEnergyRand(phonEnergy);
  G4double qp2Energy = phonEnergy - qp1Energy;
  std::pair<G4double,G4double> thePair(qp1Energy,qp2Energy);
  return thePair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Compute QP energy distribution from quasiparticle in superconductor.
// This is the same as the original KaplanQP
G4double G4CMPSCPairBreakingProcess::QPEnergyRand(G4double Energy) const
{
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy - E') + fGapEnergy*fGapEnergy)
  //           /
  //           sqrt((E'*E' - fGapEnergy*fGapEnergy) *
  //                ((Energy - E')*(Energy - E') - fGapEnergy*fGapEnergy));
  // The shape of the PDF is like a U, so the max values are at the endpoints:
  // E' = fGapEnergy and E' = Energy - fGapEnergy
  
  // Add buffer so first/last bins don't give zero denominator in pdfSum
  
  const G4double BUFF = 10000.; //REL used to be 1000
  G4double xmin = fGapEnergy + (Energy - 2. * fGapEnergy) / BUFF;
  G4double xmax = fGapEnergy + (Energy - 2. * fGapEnergy) * (BUFF - 1.) / BUFF;
  G4double ymax = QPEnergyPDF(Energy, xmin);

  G4double xtest = 0., ytest = ymax;
  do {
    ytest = G4UniformRand() * ymax;
    xtest = G4UniformRand() * (xmax - xmin) + xmin;
  } while (ytest > QPEnergyPDF(Energy, xtest));

  return xtest;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// QP energy "pdf"
// This is the same as the original KaplanQP
G4double G4CMPSCPairBreakingProcess::QPEnergyPDF(G4double E, G4double x) const
{
    const G4double gapsq = fGapEnergy * fGapEnergy;
    return ((x * (E - x) + gapsq) / sqrt((x * x - gapsq) * ((E - x) * (E - x) - gapsq)));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Take the two energies for the two pairs of QPs and actually generate secondaries
//with those energies. 
void G4CMPSCPairBreakingProcess::GenerateBogoliubovQPPair(std::pair<G4double,G4double> QPEnergies,
							  const G4Track& aTrack,
							  const G4Step& aStep)
						
{

  //Get the energies
  double energy1 = QPEnergies.first;
  double energy2 = QPEnergies.second;

  //For QPs, velocity needs to be artificially small. This is because QPs diffuse rather than propagate ballistically in most SC media, and we need
  //a special technique to handle that diffusion (and its competition with other QP processes).
  double vel1Mag = 1E-18 * CLHEP::m / CLHEP::s;
  double vel2Mag = 2E-18 * CLHEP::m / CLHEP::s; //Use this to clearly define QP ID for ID'ing which QP undergoes recombination, in case tracking IDs are weird REL
    
  //Create direction vectors for the QPs
  //For lack of something more physical, (i.e. just to test for now), we make random
  G4ThreeVector vel1 = G4RandomDirection()*vel1Mag;
  G4ThreeVector vel2 = G4RandomDirection()*vel2Mag;
  
  // Construct the secondaries and set their attributes. Hopefully we don't have to muck with the time here.
  G4Track* sec1 = G4CMP::CreateBogoliubovQP(aTrack, energy1, vel1, aTrack.GetGlobalTime(),aTrack.GetPosition());
  G4Track* sec2 = G4CMP::CreateBogoliubovQP(aTrack, energy2, vel2, aTrack.GetGlobalTime(),aTrack.GetPosition());

  //Sanity check
  if (!sec1 || !sec2) {
    G4Exception("G4CMPSCPairBreakingProcess::GenerateBogoliubovQPPair", "SCPairBreaking001",JustWarning, "Error creating secondaries");
    return;
  }

  //Finally, create the secondaries and return
  aParticleChange.SetNumberOfSecondaries(2);
  aParticleChange.AddSecondary(sec2);
  aParticleChange.AddSecondary(sec1);
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Pass-through to G4CMPVProcess class
G4double G4CMPSCPairBreakingProcess::GetMeanFreePath(const G4Track& trk, G4double prevstep, G4ForceCondition* cond)
{
  return G4CMPVProcess::GetMeanFreePath(trk,prevstep,cond);
}
