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
/// \file library/src/G4CMPVDriftProcess.cc
/// \brief Implementation of the G4CMPVDriftProcess base class
//
// $Id$

#include "G4CMPVDriftProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4CMPValleyTrackMap.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4ProcessType.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"


// Constructor and destructor

G4CMPVDriftProcess::G4CMPVDriftProcess(const G4String& processName)
  : G4VDiscreteProcess(processName, fPhonon),
    trackIVmap(G4CMPValleyTrackMap::GetInstance()), theLattice(0),
    currentTrack(0) {}

G4CMPVDriftProcess::~G4CMPVDriftProcess() {;}


// Only applies to the known charge carriers

G4bool G4CMPVDriftProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4CMPDriftElectron::Definition() ||
	  &aPD==G4CMPDriftHole::Definition() );
}


// Initialize wave vectors for currently active track(s)

void G4CMPVDriftProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions

  // Choose missing valley randomly (this should depend on momentum)
  if (!trackIVmap->Find(track)) 
    trackIVmap->SetValley(track, (int)(G4UniformRand()*4 + 1.0));

  currentTrack = track;			// Save for use by EndTracking

  // Fetch lattice for current track once, use in subsequent steps
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(track->GetVolume());
}

void G4CMPVDriftProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  trackIVmap->RemoveTrack(currentTrack);
  currentTrack = 0;
  theLattice = 0;
}


// Create new secondary track from phonon configuration

G4Track* G4CMPVDriftProcess::CreateSecondary(const G4ThreeVector& mom,
					     G4double energy, G4int iv) const {
  if (verboseLevel>1) {
    G4cout << GetProcessName() << " CreateSecondary valley " << iv
	   << " mom " << mom << " E " << energy << G4endl;
  }

  G4ParticleDefinition* pd = G4CMPDriftElectron::Definition();

  // Secondaries are created at the current track coordinates
  G4Track* sec = new G4Track(new G4DynamicParticle(pd, mom, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store wavevector in lookup table for future tracking
  trackIVmap->SetValley(sec, iv);

  if (verboseLevel>1) {
    G4cout << GetProcessName() << " secondary in valley "
	   << trackIVmap->GetValley(sec) << G4endl;
  }

  return sec;
}
