/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file validation/g4cmpValidation.cc
/// \brief Main program of the G4CMP/validation example
//
// $Id$
//

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "G4CMPPhysicsList.hh"
#include "G4CMPConfigManager.hh"
#include "ValidationActionInitialization.hh"
#include "ValidationConfigManager.hh"
#include "ValidationDetectorConstruction.hh"
#include "ValidationDetectorParameters.hh"

int main(int argc,char** argv)
{
 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

 G4VUserPhysicsList* physics = new G4CMPPhysicsList();
 physics->SetCuts();
 runManager->SetUserInitialization(physics);

 // Create configuration managers to ensure macro commands exist
 G4CMPConfigManager::Instance();
 ValidationConfigManager::Instance();

 
 // Set user action classes (different for Geant4 10.0)
 //
 runManager->SetUserInitialization(new ValidationActionInitialization);

 // Set mandatory initialization classes
 //
 ValidationDetectorConstruction* detector = new ValidationDetectorConstruction();
 runManager->SetUserInitialization(detector);

 
 // Visualization manager
 //
 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();
 
 // Get the pointer to the User Interface manager
 //
 G4UImanager* UImanager = G4UImanager::GetUIpointer();  

 if (argc==1)   // Define UI session for interactive mode
 {
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      ui->SessionStart();
      delete ui;
 }
 else           // Batch mode
 {
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
 }

 delete visManager;
 delete runManager;

 return 0;
}


