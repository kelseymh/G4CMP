#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UserSteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Tst1DetectorConstruction.hh"
#include "Tst1PhysicsList.hh"
#include "Tst1PrimaryGeneratorAction.hh"
#include "PhononStackingAction.hh"
#include "DriftingElectronStackingAction.hh"
#include "FET1.hh"
//#include "FET.hh"
//#include "FETMessenger.hh"
#include "FETMessenger1.hh"
//#include "FETUserRunAction.hh"
//#include "FETUserEventAction.hh"

int main(int argc,char** argv)
{
  

 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

 // Set mandatory initialization classes
 //
 Tst1DetectorConstruction* detector = new Tst1DetectorConstruction();
 runManager->SetUserInitialization(detector);
 //
 G4VUserPhysicsList* physics = new Tst1PhysicsList();
 physics->SetCuts();
 runManager->SetUserInitialization(physics);
    
 // Set user action classes
 //

 runManager->SetUserAction(new PhononStackingAction);
 runManager->SetUserAction(new DriftingElectronStackingAction);
 // runManager->SetUserAction(new PhononTrackingAction);

 G4VUserPrimaryGeneratorAction* gen_action = new Tst1PrimaryGeneratorAction();
 runManager->SetUserAction(gen_action);
// runManager->SetUserAction(new FETUserEventAction);
// runManager->SetUserAction(new FETUserRunAction);
 //
 FET1* fetsim = new FET1();
 //FET* fetsim = new FET(runManager);
 FETMessenger1* fetmes = new FETMessenger1(fetsim);





#ifdef G4VIS_USE
 // Visualization manager
 //
 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();
#endif
    
 // Initialize G4 kernel
 //
 runManager->Initialize();
  
 // Get the pointer to the User Interface manager
 //
 G4UImanager* UImanager = G4UImanager::GetUIpointer();  

 if (argc==1)   // Define UI session for interactive mode
 {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
 }
 else           // Batch mode
 {
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
 }

//G4cout << runManager->GetCurrentRun()->GetEventVector()->size() << G4endl;
// FET* fet = new FET(runManager->GetCurrentRun());
// delete fet;
#ifdef G4VIS_USE
 delete visManager;
#endif
 delete runManager;

 return 0;
}


