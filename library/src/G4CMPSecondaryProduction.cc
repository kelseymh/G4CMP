//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/include/G4CMPSecondaryProduction.hh
/// \brief Definition of the G4CMPSecondaryProduction process class.  This
///	class will be used to extend the existing Geant4 ionization
///	(and possibly other) processes to generate phonons and charge
///	carriers as secondaries.
//
// $Id$
//
// 20150306  Michael Kelsey

#include "G4CMPSecondaryProduction.hh"
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

// Constructor and destructor

G4CMPSecondaryProduction::G4CMPSecondaryProduction()
  : G4VContinuousProcess("G4CMPSecondaryProduction", fPhonon),
    ionizationEnergy(0.*eV), yieldPerPhonon(0.67), yieldPerChargePair(0.33),
    sigmaPerPhonon(0.05), sigmaPerChargePair(0.05) {
  SetProcessSubType(fSecondaryProduction);
}

G4CMPSecondaryProduction::~G4CMPSecondaryProduction() {;}


// Applies to all charged, non-resonance particles except the drift charges

G4bool G4CMPSecondaryProduction::IsApplicable(const G4ParticleDefinition& pd) {
  return (pd.GetPDGCharge() != 0 && !pd.IsShortLived() &&
	  &pd != G4CMPDriftElectron::Definition() &&
	  &pd != G4CMPDriftHole::Definition());
}


// Overload G4CMPProcessUtils function to fill energy parameters

void G4CMPSecondaryProduction::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);

  /***** FIXME:  MeanExcitationEnergy() below is 350 eV?!?  Not ~3 eV?
  G4Material* theMat = track->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4IonisParamMat* ionisPar = theMat->GetIonisation();

  // FIXME: We should use the parametrization to do fluctuations below
  ionizationEnergy = ionisPar->GetMeanExcitationEnergy();
  *****/
  ionizationEnergy = 3*eV;	// Hard-code this temporarily

  if (verboseLevel>1) {
    G4cout << GetProcessName() << " ionizationEnergy " << ionizationEnergy/eV
	   << " eV" << G4endl;
  }
}


// Use previously computed energy loss to generate secondaries

G4VParticleChange* 
G4CMPSecondaryProduction::AlongStepDoIt(const G4Track& track,
					const G4Step& stepData) {
  aParticleChange.Initialize(track); 

  // Only apply to tracks while they are in lattice-configured volumes
  G4VPhysicalVolume* trkPV = track.GetVolume();
  G4LatticePhysical* theLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(trkPV);
  if (!theLattice) return &aParticleChange;

  if (verboseLevel) G4cout << GetProcessName() << "::AlongStepDoIt" << G4endl;

  LoadDataForTrack(&track);		// Do on every step to change volumes
  AddPhonons(stepData);
  AddChargeCarriers(stepData);

  // NOTE:  This process does NOT change the track's momentum or energy
  return &aParticleChange;
}


// Use non-ionizing energy loss to generate phonons along path

void G4CMPSecondaryProduction::AddPhonons(const G4Step& stepData) {
  G4double Etotal = stepData.GetTotalEnergyDeposit();
  if (Etotal < ionizationEnergy) return;

  if (verboseLevel) G4cout << " AddPhonons " << Etotal/eV << " eV" << G4endl;

  G4int nsec = GenerateEnergyPositions(stepData, Etotal, yieldPerPhonon,
				       sigmaPerPhonon);

  if (verboseLevel>1) G4cout << " Adding " << nsec << " phonons" << G4endl;
  aParticleChange.SetNumberOfSecondaries(nsec);

  // Distribute generated phonons along trajectory
  G4ThreeVector dir;
  G4Track* aSec = 0;
  for (size_t i=0; i<energyPosList.size(); i++) {
    dir = G4RandomDirection();
    aSec = CreatePhonon(ChoosePolarization(), dir, energyPosList[i].first,
			energyPosList[i].second);

    if (verboseLevel>2) {
      G4cout << " Created secondary " << i << " "
	     << aSec->GetParticleDefinition()->GetParticleName()
	     << "\n track along " << aSec->GetMomentumDirection()
	     << " (length " << aSec->GetMomentumDirection().mag() << ")"
	     << "\n E " << aSec->GetKineticEnergy()/eV << " eV,"
	     << "\n @ " << aSec->GetPosition() << G4endl;
    }

    aParticleChange.AddSecondary(aSec);
  }
}


// Use ionizing energy loss to generate charge pairs along path

void G4CMPSecondaryProduction::AddChargeCarriers(const G4Step& stepData) {
  G4double Etotal = stepData.GetTotalEnergyDeposit();
  if (Etotal < ionizationEnergy) return;

  if (verboseLevel)
    G4cout << " AddChargeCarriers " << Etotal/eV << " eV" << G4endl;

  G4int nsec = GenerateEnergyPositions(stepData, Etotal, yieldPerChargePair,
				       sigmaPerChargePair);

  if (verboseLevel>1) G4cout << " Adding " << nsec << " e/h pairs" << G4endl;
  aParticleChange.SetNumberOfSecondaries(2*nsec);

  // Distribute generated phonons along trajectory
  G4Track* aSec = 0;
  G4ThreeVector dir;
  for (size_t i=0; i<energyPosList.size(); i++) {
    dir = G4RandomDirection();		// Produce back-to-back

    // FIXME:  Energy should be split to balance momenta!
    aSec = CreateChargeCarrier(-1, ChooseValley(), energyPosList[i].first/2.,
			       dir, energyPosList[i].second);
    if (verboseLevel>2) {
      G4cout << " Created secondary " << 2*i << " "
	     << aSec->GetParticleDefinition()->GetParticleName()
	     << "\n along " << aSec->GetMomentumDirection()
	     << " (length " << aSec->GetMomentumDirection().mag() << ")"
	     << "\n E " << aSec->GetKineticEnergy()/eV << " eV,"
	     << "\n @ " << aSec->GetPosition() << G4endl;
    }

    aParticleChange.AddSecondary(aSec);

    aSec = CreateChargeCarrier(1, 0, energyPosList[i].first/2.,
			       -dir, energyPosList[i].second);
    if (verboseLevel>2) {
      G4cout << " Created secondary " << 2*i+1 << " "
	     << aSec->GetParticleDefinition()->GetParticleName()
	     << "\n track along " << aSec->GetMomentumDirection()
	     << " (length " << aSec->GetMomentumDirection().mag() << ")"
	     << "\n E " << aSec->GetKineticEnergy()/eV << " eV,"
	     << "\n @ " << aSec->GetPosition() << G4endl;
    }

    aParticleChange.AddSecondary(aSec);
  }
}


// Generate energy steps along step trajectory (straight line only!)

G4int G4CMPSecondaryProduction::GenerateEnergyPositions(const G4Step& stepData,
							G4double Etotal,
							G4double yield, 
							G4double sigma) {
  if (verboseLevel>1) {
    G4cout << " GenerateEnergyPositions Etotal " << Etotal/eV << " eV"
	   << " yield " << yield << " sigma " << sigma << G4endl;
  }

  G4int nsec = Etotal / ionizationEnergy;	// Average number of secondaries

  // Get average distance between secondaries along (straight) trajectory
  G4ThreeVector prePos  = stepData.GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = stepData.GetPostStepPoint()->GetPosition();
  G4ThreeVector traj = postPos - prePos;

  G4double length = (postPos-prePos).mag();
  G4double dr = length / G4double(nsec);

  if (verboseLevel>1) {
    G4cout << " About " << nsec << " secondaries along " << length/mm << " mm"
	   << G4endl;
  }

  // Minimum energy and spacing values to stop iteration
  G4double emin = ionizationEnergy * (1. - 2.*sigma);
  G4double lmin = 0.5 * dr;

  // Build list of energy-position pairs for secondaries
  energyPosList.clear();
  energyPosList.reserve(nsec);		// Assume we get the average

  EPosPair theEPos(0.,prePos);		// Buffers to fill secondary steps
  G4double ystep=0., substep=0.;
  do {
    ystep = G4RandGauss::shoot(yield, sigma);
    substep = G4RandExponential::shoot(dr);	// Distance to next step

    theEPos.first = ionizationEnergy * ystep;
    theEPos.second += substep*traj;
    energyPosList.push_back(theEPos);

    Etotal -= theEPos.first/yield;		// May end up negative at end
    length -= substep;
  } while (Etotal > emin && length > lmin);

  // Adjust energies and lengths to absorb residuals
  nsec = (G4int)energyPosList.size();
  G4double eadj = Etotal / G4double(nsec);
  G4double ladj = length / G4double(nsec);

  for (G4int i=0; i<nsec; i++) {
    energyPosList[i].first  += eadj;
    energyPosList[i].second += ladj*traj;
  }

  // Resulting list returned for use in building secondaries
  return nsec;
}


// Calculate step limit for Along Step (not needed here)

G4double 
G4CMPSecondaryProduction::GetContinuousStepLimit(const G4Track& aTrack,
						 G4double  previousStepSize,
						 G4double  currentMinimumStep,
						 G4double& currentSafety) {
  return DBL_MAX;	// This should prevent step-limiting here
}

