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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1AnalysisConfig.hh"
// #include "B1Run.hh"

#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1100
#include "G4AnalysisManager.hh"
#else // G4VERSION_NUMBER >= 1100
#include "g4root.hh"
// #include "g4csv.hh"
#endif // G4VERSION_NUMBER >= 1100

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Ntuple 0: interface hits
  analysisManager->CreateNtuple("InterfaceHits", "Particles reaching Si/Al interface");
  analysisManager->CreateNtupleIColumn("eventID");      // 0
  analysisManager->CreateNtupleIColumn("trackID");      // 1
  analysisManager->CreateNtupleIColumn("speciesID");    // 2
  analysisManager->CreateNtupleDColumn("weight");       // 3
  analysisManager->CreateNtupleDColumn("energy_eV");    // 4
  analysisManager->CreateNtupleDColumn("vperp_mm_ns");  // 5
  analysisManager->CreateNtupleDColumn("time_ns");      // 6
  analysisManager->CreateNtupleDColumn("x_mm");         // 7
  analysisManager->CreateNtupleDColumn("y_mm");         // 8
  analysisManager->CreateNtupleDColumn("z_mm");         // 9
  analysisManager->FinishNtuple(0);

  // Ntuple 1: electron event summary
  analysisManager->CreateNtuple("EventSummary", "Per-event electron summary");
  analysisManager->CreateNtupleIColumn("eventID");         // col 0
  analysisManager->CreateNtupleDColumn("nInSi");           // col 1
  analysisManager->CreateNtupleDColumn("nReached");        // col 2
  analysisManager->CreateNtupleDColumn("nTrapped");        // col 3
  analysisManager->CreateNtupleDColumn("reachFraction");   // col 4
  analysisManager->FinishNtuple(1);

  // Ntuple 2: phonon event summary
  analysisManager->CreateNtuple("PhononEventSummary", "Per-event phonon interface summary");
  analysisManager->CreateNtupleIColumn("eventID");                       // col 0
  analysisManager->CreateNtupleDColumn("nPhononsInterfaceWeighted");     // col 1
  analysisManager->CreateNtupleDColumn("ePhononsInterfaceWeighted_eV");  // col 2
  analysisManager->CreateNtupleDColumn("nPhononsAbove2DeltaWeighted");   // col 3
  analysisManager->CreateNtupleDColumn("ePhononsAbove2DeltaWeighted_eV");// col 4
  analysisManager->FinishNtuple(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{
  fEdep = 0.;

  auto analysisManager = G4AnalysisManager::Instance();
  auto config = B1AnalysisConfig::Instance();

  G4cout << "Opening analysis output file: "
         << config->GetOutputFileName() << G4endl;

  analysisManager->OpenFile(config->GetOutputFileName());
}

void B1RunAction::EndOfRunAction(const G4Run*)
{
  G4cout << G4endl
         << "Total Edep eV " << fEdep / eV << G4endl;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

void B1RunAction::AddEdep(G4double edep)
{
  fEdep += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

