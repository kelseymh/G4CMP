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
/// \file library/include/G4CMPIonisationWrapper.hh
/// \brief Definition of the G4CMPIonisationWrapper process class.  This
///	class will be used to extend the existing Geant4 ionization
///	(and possibly other) processes to generate phonons and charge
///	carriers as secondaries.
//
// $Id$
//
// 20150306  Michael Kelsey

#include "G4CMPIonisationWrapper.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"

// Constructor and destructor

G4CMPIonisationWrapper::G4CMPIonisationWrapper(G4VProcess* ionizProc)
  : G4WrapperProcess(ionizProc->GetProcessName(), ionizProc->GetProcessType()),
    energyPerPhonon(1e-3*eV), sigmaEPerPhonon(1e-4*eV),
    energyPerChargePair(1.*eV), sigmaEPerChargePair(0.1*eV) {
  RegisterProcess(ionizProc);
}

G4CMPIonisationWrapper::~G4CMPIonisationWrapper() {;}


// Actions:  Let process do its work, use energy loss to generate secondaries

G4VParticleChange* 
G4CMPIonisationWrapper::PostStepDoIt(const G4Track& track,
			             const G4Step&  stepData) {
  LoadDataForTrack(&track);	// Get lattice information

  G4VParticleChange* theChange = pRegProcess->PostStepDoIt(track, stepData);
  AddPhonons(theChange, stepData);
  AddChargeCarriers(theChange, stepData);

  ReleaseTrack();

  return theChange;
}

G4VParticleChange* 
G4CMPIonisationWrapper::AlongStepDoIt(const G4Track& track,
				      const G4Step& stepData) {
  LoadDataForTrack(&track);	// Get lattice information

  G4VParticleChange* theChange = pRegProcess->AlongStepDoIt(track, stepData);
  AddPhonons(theChange, stepData);
  AddChargeCarriers(theChange, stepData);

  ReleaseTrack();

  return theChange;
}


// Use non-ionizing energy loss to generate phonons along path

void G4CMPIonisationWrapper::AddPhonons(G4VParticleChange* theChange,
					const G4Step& stepData) {
  G4double Etotal = theChange->GetNonIonizingEnergyDeposit();
  if (Etotal < energyPerPhonon) return;

  if (verboseLevel) G4cout << "G4CMPIonisationWrapper::AddPhonons" << G4endl;

  GenerateEnergyPositions(stepData, Etotal, energyPerPhonon, sigmaEPerPhonon);

  // Distribute generated phonons along trajectory
  G4ThreeVector prePos  = stepData.GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = stepData.GetPostStepPoint()->GetPosition();
  G4ThreeVector traj = postPos - prePos;

  G4ThreeVector pos = prePos;
  G4Track* aSec = 0;
  for (size_t i=0; i<energyPosList.size(); i++) {
    aSec = CreatePhonon(ChoosePolarization(), G4RandomDirection(),
			energyPosList[i].first, pos);
    theChange->AddSecondary(aSec);

    pos += energyPosList[i].second * traj;
  }
}


// Use ionizing energy loss to generate charge pairs along path

void G4CMPIonisationWrapper::AddChargeCarriers(G4VParticleChange* theChange,
					       const G4Step& stepData) {
  G4double Etotal = (theChange->GetLocalEnergyDeposit()
		     - theChange->GetNonIonizingEnergyDeposit());
  if (Etotal < energyPerChargePair) return;

  if (verboseLevel)
    G4cout << "G4CMPIonisationWrapper::AddChargeCarriers" << G4endl;

  GenerateEnergyPositions(stepData, Etotal, energyPerChargePair,
			  sigmaEPerChargePair);


  // Distribute generated phonons along trajectory
  G4ThreeVector prePos  = stepData.GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = stepData.GetPostStepPoint()->GetPosition();
  G4ThreeVector traj = postPos - prePos;

  G4ThreeVector pos = prePos;
  G4ThreeVector dir;
  G4Track* aSec = 0;
  for (size_t i=0; i<energyPosList.size(); i++) {
    dir = G4RandomDirection();		// Produce with back-to-back momentum
    aSec = CreateChargeCarrier(-1, ChooseValley(), energyPosList[i].first/2.,
			       dir, pos);
    theChange->AddSecondary(aSec);

    aSec = CreateChargeCarrier(1, ChooseValley(), energyPosList[i].first/2.,
			       -dir, pos);
    theChange->AddSecondary(aSec);

    pos += energyPosList[i].second * traj;
  }
}


// Generate energy steps along step trajectory (straight line only!)

void G4CMPIonisationWrapper::GenerateEnergyPositions(const G4Step& stepData,
						     G4double Etotal,
						     G4double Eunit, 
						     G4double sigmaE) {
  if (verboseLevel>1) {
    G4cout << " GenerateEnergyPositions Etotal " << Etotal/eV
	   << " Eunit " << Eunit/eV << " sigmaE " << sigmaE/eV
	   << " [eV]" << G4endl;
  }

  G4int nsec = Etotal / Eunit;		// Average number of secondaries

  // Get average distance between secondaries along (straight) trajectory
  G4ThreeVector prePos  = stepData.GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = stepData.GetPostStepPoint()->GetPosition();
  G4double length = (postPos-prePos).mag();
  G4double dr = length / G4double(nsec);

  // Minimum energy and spacing values to stop iteration
  G4double emin = Eunit - 2.*sigmaE;
  G4double lmin = 0.5 * dr;

  // Build list of energy-position pairs for secondaries
  EPosPair theEPos(0.,0.);
  energyPosList.clear();
  energyPosList.reserve(nsec);		// Assume we get the average
  do {
    theEPos.first = G4RandGauss::shoot(Eunit, sigmaE);
    theEPos.second = G4RandExponential::shoot(dr);
    
    Etotal -= theEPos.first;		// May end up negative at end
    length -= theEPos.second;
  } while (Etotal > emin && length > lmin);

  // Adjust energies and lengths to absorb residuals
  nsec = (G4int)energyPosList.size();
  G4double eadj = Etotal / G4double(nsec);
  G4double ladj = length / G4double(nsec);

  for (G4int i=0; i<nsec; i++) {
    energyPosList[i].first += eadj;
    energyPosList[i].second += ladj;
  }

  // Resulting list returned for use in building secondaries
}
