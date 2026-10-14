/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file PhononSensitivity.cc
/// \brief Implementation of the PhononSensitivity class
///
/// Definition of hits, hits collection, output scheme

#include "PhononSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "PhononConfigManager.hh"
#include <fstream>


PhononSensitivity::PhononSensitivity(G4String name) :
  G4CMPElectrodeSensitivity(name), fileName("") {
  SetOutputFile(PhononConfigManager::GetHitOutput());
}

PhononSensitivity::~PhononSensitivity() {
  if (output.is_open()) output.close();
  if (!output.good()) {
    G4cerr << "Error closing output file, " << fileName << ".\n"
           << "Expect bad things like loss of data.";
  }
}

void PhononSensitivity::EndOfEvent(G4HCofThisEvent* HCE) {
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  auto* hitCol = static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(HCID));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();

  G4RunManager* runMan = G4RunManager::GetRunManager();

  if (output.good()) {
    for (G4CMPElectrodeHit* hit : *hitVec) {
      output << runMan->GetCurrentRun()->GetRunID() << ','
             << runMan->GetCurrentEvent()->GetEventID() << ','
             << hit->GetTrackID() << ','
             << hit->GetParticleName() << ','
             << hit->GetStartEnergy()/eV << ','
             << hit->GetStartPosition().getX()/m << ','
             << hit->GetStartPosition().getY()/m << ','
             << hit->GetStartPosition().getZ()/m << ','
             << hit->GetStartTime()/ns << ','
             << hit->GetEnergyDeposit()/eV << ','
             << hit->GetWeight() << ','
             << hit->GetFinalPosition().getX()/m << ','
             << hit->GetFinalPosition().getY()/m << ','
             << hit->GetFinalPosition().getZ()/m << ','
             << hit->GetFinalTime()/ns << '\n';
    }
  }
}

void PhononSensitivity::SetOutputFile(const G4String &fn) {
  if (fileName != fn) {
    if (output.is_open()) output.close();
    fileName = fn;
    output.open(fileName, std::ios_base::app);
    if (!output.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening output file " << fileName;
      G4Exception("PhononSensitivity::SetOutputFile", "PhonSense003",
                  FatalException, msg);
      output.close();
    } else {
      output << "Run_ID,Event_ID,Track_ID,Particle_Name,Start_Energy_[eV],"
             << "Start_X_[m],Start_Y_[m],Start_Z_[m],Start_Time_[ns],"
             << "Energy_Deposited_[eV],Track_Weight,End_X_[m],End_Y_[m],End_Z_[m],"
             << "Final_Time_[ns]\n";
    }
  }
}

G4bool PhononSensitivity::IsHit(const G4Step* step,
                                const G4TouchableHistory*) const {
  /* Phonons tracks are sometimes killed at the boundary in order to spawn new
   * phonon tracks. These tracks that are killed deposit no energy and should
   * not be picked up as hits.
   */
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition* particle = track->GetDefinition();

  G4bool correctParticle = particle == G4PhononLong::Definition() ||
                           particle == G4PhononTransFast::Definition() ||
                           particle == G4PhononTransSlow::Definition();

  G4bool correctStatus = step->GetTrack()->GetTrackStatus() == fStopAndKill &&
                         postStepPoint->GetStepStatus() == fGeomBoundary &&
                         step->GetNonIonizingEnergyDeposit() > 0.;

  return correctParticle && correctStatus;
}
