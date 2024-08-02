/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

//20240110 Israel Hernandez -- Illinois Institute of Technology
#include "Caustic_PhononSensitivity.hh"
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
#include "Caustic_PhononConfigManager.hh"
#include <fstream>
//20240110 Israel Hernandez -- Illinois Institute of Technology

Caustic_PhononSensitivity::Caustic_PhononSensitivity(G4String name) :
  G4CMPElectrodeSensitivity(name), fileName("") {
  SetOutputFile(Caustic_PhononConfigManager::GetHitOutput());
}



Caustic_PhononSensitivity::~Caustic_PhononSensitivity() {
  if (output.is_open()) output.close();
  if (!output.good()) {
    G4cerr << "Error closing output file, " << fileName << ".\n"
           << "Expect bad things like loss of data.";
  }
}

void Caustic_PhononSensitivity::EndOfEvent(G4HCofThisEvent* HCE) {
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  auto* hitCol = static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(HCID));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();

  G4RunManager* runMan = G4RunManager::GetRunManager();

  if (output.good()) {
    //Saving in a txt file the Final Phonon Position.
    for (G4CMPElectrodeHit* hit : *hitVec) {
        output << hit->GetParticleName() << '\t'
        << hit->GetFinalPosition().getX()/m << '\t'
        << hit->GetFinalPosition().getY()/m << '\t'
        << hit->GetFinalPosition().getZ()/m << '\n';

    }
  }
}

void Caustic_PhononSensitivity::SetOutputFile(const G4String &fn) {
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

    }
  }
}

G4bool Caustic_PhononSensitivity::IsHit(const G4Step* step,
                                const G4TouchableHistory*) const {
  /* Phonons tracks are sometimes killed at the boundary in order to spawn new
   * phonon tracks. These tracks that are killed deposit no energy and should
   * not be picked up as hits.
   */
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition* particle = track->GetDefinition();
// I include this only to save the data on  the Aluminum Detector
  const G4TouchableHandle touch1 = postStepPoint->GetTouchableHandle();
  const G4VPhysicalVolume* volume = touch1->GetVolume();
  const G4String name = volume->GetName();

  G4bool correctParticle = particle == G4PhononLong::Definition() ||
                           particle == G4PhononTransFast::Definition() ||
                           particle == G4PhononTransSlow::Definition();
                           // I added the additional condition to save only the phonon that impacts  the Aluminum Sensor
  G4bool correctStatus = step->GetTrack()->GetTrackStatus() == fStopAndKill &&
                         postStepPoint->GetStepStatus() == fGeomBoundary &&
                         step->GetNonIonizingEnergyDeposit() > 0. && name=="fBolometerPhysical";

  return correctParticle && correctStatus;
}
