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
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Version.hh"

// Start new for application
#if G4VERSION_NUMBER >= 1100
#include "G4AnalysisManager.hh"
#else // G4VERSION_NUMBER >= 1100
#include "g4root.hh"
// #include "g4csv.hh"
#endif // G4VERSION_NUMBER >= 1100
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4CMPDriftElectron.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "B1AnalysisConfig.hh"
// End new for application

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
// Start new for application
  fSiVolume(nullptr),
  fAlVolume(nullptr)
// End new for application
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Begin new SteppingAction for application
void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fSiVolume || !fAlVolume) {
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    fScoringVolume = detectorConstruction->GetScoringVolume();
    fSiVolume      = detectorConstruction->GetSiLogicalVolume();
    fAlVolume      = detectorConstruction->GetAlLogicalVolume();
  }

  G4Track* track = step->GetTrack();
  if (!track) return;

  G4StepPoint* prePoint  = step->GetPreStepPoint();
  G4StepPoint* postPoint = step->GetPostStepPoint();
  if (!prePoint || !postPoint) return;

  G4VPhysicalVolume* prePV  = prePoint->GetPhysicalVolume();
  G4VPhysicalVolume* postPV = postPoint->GetPhysicalVolume();

  G4LogicalVolume* preLV  = prePV  ? prePV->GetLogicalVolume()  : nullptr;
  /*** Avoid compiler warning about "unused variable"
  G4LogicalVolume* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;
  ***/

  G4String preName  = prePV  ? prePV->GetName()  : "";
  G4String postName = postPV ? postPV->GetName() : "";

  G4int trackID = track->GetTrackID();
                                  
  //Fetch weight associated with track from downsampling
  G4double weight = track->GetWeight();

  // Restrict this tutorial stage to electrons only
  G4String particleName = track->GetDefinition()->GetParticleName();
  G4bool acceptTrack = false;
  G4int speciesID = -1;

  auto config = B1AnalysisConfig::Instance();

  if (config->IsElectronMode()) {
    if (particleName == "G4CMPDriftElectron" || particleName == "G4CMPDriftHole") {
      acceptTrack = true;
      speciesID = 0;
    }
  }
  else if (config->IsPhononMode()) {
    if (particleName == "phononL") {
      acceptTrack = true;
      speciesID = 1;
    }
    else if (particleName == "phononTF") {
      acceptTrack = true;
      speciesID = 2;
    }
    else if (particleName == "phononTS") {
      acceptTrack = true;
      speciesID = 3;
    }
  }

  if (!acceptTrack) return;

  // ---------------------------------------------------------
  // 1. Count electrons that exist in Si
  //    Use first step in Si as "created/entered Si" bookkeeping.
  // ---------------------------------------------------------

  if (config->IsElectronMode()) {
    if (preName == "Si_box" || postName == "Si_box") {
      fEventAction->CountElectronInSi(trackID, weight);
    }
  }
  // ---------------------------------------------------------
  // 2. Detect crossing from Si to Al
  // ---------------------------------------------------------

  if (preName == "Si_box" && postName == "Al_box") {
    G4double energy   = prePoint->GetKineticEnergy();
    G4double time     = prePoint->GetGlobalTime();
    G4ThreeVector pos = prePoint->GetPosition();

    G4ThreeVector momDir = prePoint->GetMomentumDirection();
    G4double speed = track->GetVelocity();
    G4ThreeVector interfaceNormal(0., 0., 1.);
    G4double vperp = speed * momDir.dot(interfaceNormal);

    fEventAction->CountElectronReachedInterface(trackID, weight);
    fEventAction->AddInterfaceElectronKE(energy);
    fEventAction->AddInterfaceElectronVperp(vperp);
    fEventAction->AddInterfaceElectronTime(time);
    fEventAction->AddInterfaceElectronPosition(pos);

    auto analysisManager = G4AnalysisManager::Instance();
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleIColumn(0, 1, trackID);
    analysisManager->FillNtupleIColumn(0, 2, speciesID);
    analysisManager->FillNtupleDColumn(0, 3, weight);
    analysisManager->FillNtupleDColumn(0, 4, energy / eV);
    analysisManager->FillNtupleDColumn(0, 5, vperp / (mm/ns));
    analysisManager->FillNtupleDColumn(0, 6, time / ns);
    analysisManager->FillNtupleDColumn(0, 7, pos.x() / mm);
    analysisManager->FillNtupleDColumn(0, 8, pos.y() / mm);
    analysisManager->FillNtupleDColumn(0, 9, pos.z() / mm);
    analysisManager->AddNtupleRow(0);

    // If we are in phonon mode, accumulate phonon event-level interface quantities
    if (config->IsPhononMode()) {
      fEventAction->AddPhononInterfaceHit(weight, energy);
  }
  }
  // ---------------------------------------------------------
  // 3. Detect electron termination in Si
  // ---------------------------------------------------------
  if (config->IsElectronMode()) {
    G4TrackStatus status = track->GetTrackStatus();
    if (postName == "Si_box" &&
        (status == fStopAndKill || status == fKillTrackAndSecondaries)) {
      fEventAction->CountElectronTrappedInSi(trackID, weight);
    }
  }

  // ---------------------------------------------------------
  // 4. Retain old Al energy deposition scoring 
  // ---------------------------------------------------------
  if (preLV == fScoringVolume) {
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  }
}
// End new SteppingAction for application


