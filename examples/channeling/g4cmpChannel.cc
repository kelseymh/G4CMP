/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4RunManager.hh"

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "ChannelingDetectorConstruction.hh"
#include "ChannelingPrimaryGeneratorAction.hh"
#include "ChannelingPhysicsList.hh"
#include "ChannelingTrackingAction.hh"
#include "ChannelingStackingAction.hh"
#include "ChannelingEventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;

    // Activate UI-command base scorer
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);
    
    // Choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    
    // Set mandatory initialization classes
    G4VUserDetectorConstruction* detector = new ChannelingDetectorConstruction;
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new ChannelingPhysicsList());

    // Set user action classes
    runManager->SetUserAction(new ChannelingPrimaryGeneratorAction());
    runManager->SetUserAction(new ChannelingEventAction());
    runManager->SetUserAction(new ChannelingStackingAction());
    runManager->SetUserAction(new ChannelingTrackingAction());
       
    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();  
    
    if(argc!=1) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    else {
        // Define visualization and UI terminal for interactive mode
#ifdef G4VIS_USE
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
#endif
        
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
        ui->SessionStart();
        delete ui;
#endif
        
#ifdef G4VIS_USE
        delete visManager;
#endif
    }
    
    // Job termination
    delete runManager;
    
    return 0;
}
