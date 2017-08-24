/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file charge/g4cmpIonize.cc
/// \brief Alternate G4CMP/charge example, demonstrating phonon, e/h production
//
// $Id$
//
// 20170816  Add example-specific configuration manager

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "ChargeActionInitialization.hh"
#include "ChargeDetectorConstruction.hh"
#include "ChargeConfigManager.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPPhysics.hh"
#include "G4CMPStackingAction.hh"
#include "FTFP_BERT.hh"


int main(int argc,char** argv) {
  // Construct the run manager
  //
  G4RunManager* runManager = new G4RunManager;
  
  // Set mandatory initialization classes
  //
  ChargeDetectorConstruction* detector = new ChargeDetectorConstruction;
  runManager->SetUserInitialization(detector);
  
  FTFP_BERT* physics = new FTFP_BERT;
  physics->RegisterPhysics(new G4CMPPhysics);
  physics->SetCuts();
  runManager->SetUserInitialization(physics);
  
  // Set user action classes (different for Geant4 10.0)
  //
  runManager->SetUserInitialization(new ChargeActionInitialization);
  
  // Create G4CMP configuration manager to ensure macro commands exist
  G4CMPConfigManager::Instance();
  ChargeConfigManager::Instance();
  
#ifdef G4VIS_USE
  // Visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();
#endif
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  
  if (argc==1) {   // Define UI session for interactive mode
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    ui->SessionStart();
    delete ui;
#endif
  } else {           // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  
  return 0;
}
