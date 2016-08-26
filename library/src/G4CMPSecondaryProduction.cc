/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSecondaryProduction.hh
/// \brief Definition of the G4CMPSecondaryProduction process class.  This
///	class will be used to extend the existing Geant4 ionization
///	(and possibly other) processes to generate phonons and charge
///	carriers as secondaries.
//
// $Id$
//
// 20150306  Michael Kelsey
// 20160825  Replace implementation with use of G4CMPEnergyPartition

#include "G4CMPSecondaryProduction.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPProcessSubType.hh"
#include "G4IonisParamMat.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include <algorithm>
#include <vector>


// Constructor and destructor

G4CMPSecondaryProduction::G4CMPSecondaryProduction()
  : G4VContinuousProcess("G4CMPSecondaryProduction", fPhonon),
    partitioner(new G4CMPEnergyPartition) {
  SetProcessSubType(fSecondaryProduction);
}

G4CMPSecondaryProduction::~G4CMPSecondaryProduction() {
  delete partitioner;
}


// Applies to all charged, non-resonance particles except the drift charges

G4bool G4CMPSecondaryProduction::IsApplicable(const G4ParticleDefinition& pd) {
  return (pd.GetPDGCharge() != 0 && !pd.IsShortLived() &&
	  &pd != G4CMPDriftElectron::Definition() &&
	  &pd != G4CMPDriftHole::Definition());
}


// Overload G4CMPProcessUtils function to fill energy parameters

void G4CMPSecondaryProduction::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);
  *partitioner = *this;		// Copies lattice and other important info

  partitioner->SetVerboseLevel(verboseLevel);
}


// Use previously computed energy loss to generate secondaries

G4VParticleChange* 
G4CMPSecondaryProduction::AlongStepDoIt(const G4Track& track,
					const G4Step& stepData) {
  aParticleChange.Initialize(track); 

  // Only apply to tracks while they are in lattice-configured volumes
  G4VPhysicalVolume* trkPV = track.GetVolume();
  G4LatticePhysical* lattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(trkPV);
  if (!lattice) return &aParticleChange;

  if (verboseLevel) G4cout << GetProcessName() << "::AlongStepDoIt" << G4endl;

  LoadDataForTrack(&track);		// Do on every step to change volumes
  AddSecondaries(stepData);

  // NOTE:  This process does NOT change the track's momentum or energy
  return &aParticleChange;
}


// Use energy loss to generate phonons and charge carriers along path

void G4CMPSecondaryProduction::AddSecondaries(const G4Step& stepData) {
  G4double eTotal = stepData.GetTotalEnergyDeposit();
  G4double eNIEL  = stepData.GetNonIonizingEnergyDeposit();

  if (verboseLevel)
    G4cout << " AddSecondaries " << eTotal/eV << " eV" << G4endl;

  // Configure energy partitioning for EM, nuclear, or pre-determined energy
  if (eNIEL <= 0.) {
    if (GetCurrentParticle()->GetParticleType() == "nucleus")
      partitioner->NuclearRecoil(eTotal);
    else
      partitioner->Ionization(eTotal);
  } else {
    partitioner->DoPartition(eTotal-eNIEL, eTotal);
  }

  partitioner->GetSeconaries(theSecs);
  size_t nsec = theSecs.size();

  if (verboseLevel>1) G4cout << " Adding " << nsec << " secondaries" << G4endl;
  aParticleChange.SetNumberOfSecondaries(nsec);

  // Distribute generated particles along (straight line) trajectory
  std::random_shuffle(theSecs.begin(), theSecs.end(), RandomIndex);

  GeneratePositions(stepData, nsec);
  for (size_t i=0; i<theSecs.size(); i++) {
    theSecs[i].SetPosition(posSecs[i]);
    aParticleChange.AddSecondary(theSecs[i]);

    if (verboseLevel>2) {
      G4Track* aSec = theSecs[i];
      G4cout << " Created secondary " << i << " "
	     << aSec->GetParticleDefinition()->GetParticleName()
	     << "\n track along " << aSec->GetMomentumDirection()
	     << " (length " << aSec->GetMomentumDirection().mag() << ")"
	     << "\n E " << aSec->GetKineticEnergy()/eV << " eV,"
	     << "\n @ " << aSec->GetPosition() << G4endl;
    }
  }
}


// Generate intermediate points along step trajectory (straight line only!)
// NOTE:  For MSC type deposition, these points ought to be a random walk

void G4CMPSecondaryProduction::GeneratePositions(const G4Step& stepData,
						 size_t nsec) {
  if (verboseLevel>1) G4cout << " GeneratePositions " << nsec << G4endl;

  // Get average distance between secondaries along (straight) trajectory
  G4ThreeVector prePos  = stepData.GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = stepData.GetPostStepPoint()->GetPosition();
  G4ThreeVector traj = postPos - prePos;
  G4ThreeVector tdir = traj.unit();

  G4double length = (postPos-prePos).mag();
  G4double dl = length / G4double(nsec);
  G4double sigl = dr/6.;

  if (verboseLevel>1)
    G4cout << " Choosing positions along " << length/mm << " mm" << G4endl;

  posSecs.clear();
  posSecs.reserve(nsec);

  G4double substep = 0.;
  G4ThreeVector lastPos = prePos;
  for (size_t i=0; i<nsec; i++) {
    substep = G4RandGauss::shoot(dl, sigl);
    lastPos += substep*tdir;
    posSecs.push_back(lastPos);
  }
}


// Calculate step limit for Along Step (not needed here)

G4double 
G4CMPSecondaryProduction::GetContinuousStepLimit(const G4Track& /*aTrack*/,
						 G4double  /*previousStepSize*/,
						 G4double  /*currentMinimumStep*/,
						 G4double& /*currentSafety*/) {
  return DBL_MAX;	// This should prevent step-limiting here
}


// Generate random index for shuffling secondaries

size_t G4CMPSecondaryProduction::RandomIndex(size_t n) {
  return (size_t)(n*G4UniformRand());
}
