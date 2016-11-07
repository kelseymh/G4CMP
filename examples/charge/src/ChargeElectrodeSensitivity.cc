/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeElectrodeSensitivity.hh"
#include "ChargeFETDigitizerModule.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"

#include <fstream>

ChargeElectrodeSensitivity::ChargeElectrodeSensitivity(G4String name) :
  G4VSensitiveDetector(name),
  FET(new ChargeFETDigitizerModule("FETSim")),
  hitsCollection(nullptr), fileName("") {
  SetOutputFile(G4CMPConfigManager::GetHitOutput());
  collectionName.insert("G4CMPElectrodeHit");
}

ChargeElectrodeSensitivity::ChargeElectrodeSensitivity(ChargeElectrodeSensitivity&& in) :
  G4VSensitiveDetector(std::move(in)), FET(std::move(in.FET)),
  hitsCollection(in.hitsCollection), output(std::move(in.output)),
  fileName(in.fileName) {
}

ChargeElectrodeSensitivity& ChargeElectrodeSensitivity::operator=(ChargeElectrodeSensitivity&& in) {
  // Base class members
  SensitiveDetectorName = std::move(in.SensitiveDetectorName);
  thePathName = std::move(in.thePathName);
  fullPathName = std::move(in.fullPathName);
  verboseLevel = in.verboseLevel;
  active = in.active;
  std::swap(ROgeometry, in.ROgeometry);
  std::swap(filter, in.filter);

  // Our members
  FET = std::move(in.FET);
  hitsCollection = in.hitsCollection;
  output.close();
  output = std::move(in.output);
  fileName = in.fileName;

  return *this;
}

ChargeElectrodeSensitivity::~ChargeElectrodeSensitivity() {
  if (output.is_open()) output.close();
  if (!output.good()) {
    G4cerr << "Error closing output file, " << fileName << ".\n"
           << "Expect bad things like loss of data.";
  }
}

void ChargeElectrodeSensitivity::Initialize(G4HCofThisEvent* HCE) {
  hitsCollection = new G4CMPElectrodeHitsCollection(SensitiveDetectorName,
                                                    collectionName[0]);
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  HCE->AddHitsCollection(HCID, hitsCollection);
}

void ChargeElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* HCE) {
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
  FET->Digitize();
}

void ChargeElectrodeSensitivity::SetOutputFile(const G4String &fn) {
  if (fileName != fn) {
    if (output.is_open()) output.close();
    fileName = fn;
    output.open(fileName, std::ios_base::app);
    if (!output.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening output file, " << fileName << ".\n"
          << "Will continue simulation.";
      G4Exception("ChargeElectrodeSensitivity::SetOutputFile", "Charge003",
                  JustWarning, msg);
      output.close();
    } else {
      output << "Run ID,Event ID,Track ID,Particle Name,Start Energy [eV],"
             << "Start X [m],Start Y [m],Start Z [m],Start Time [ns],"
             << "Energy Deposited [eV],Track Weight,End X [m],End Y [m],End Z [m],"
             << "Final Time [ns]\n";
    }
  }
}

G4bool ChargeElectrodeSensitivity::ProcessHits(G4Step* aStep,
                                               G4TouchableHistory* ROhist) {
  if (IsHit(aStep, ROhist)) {
    auto hit = new G4CMPElectrodeHit;
    G4CMP::FillHit(aStep, hit); // Mutates hit
    hitsCollection->insert(hit);
  }

  return true;
}

G4bool ChargeElectrodeSensitivity::IsHit(const G4Step* step,
                                        const G4TouchableHistory*) const {
  // Charge carriers do not deposit energy when they land on an electrode.
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition* particle = track->GetDefinition();

  G4bool correctParticle = particle == G4CMPDriftElectron::Definition() ||
                           particle == G4CMPDriftHole::Definition();

  G4bool correctStatus = step->GetTrack()->GetTrackStatus() == fStopAndKill;
  //G4bool correctStatus = step->GetTrack()->GetTrackStatus() == fStopAndKill &&
  //                       postStepPoint->GetStepStatus() == fGeomBoundary;

  return correctParticle && correctStatus;
}
