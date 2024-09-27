/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file phonon/g4cmpPhonon.cc
/// \brief Main program of the G4CMP/phonon example
//
// $Id$
//
// 20140509  Add conditional code for Geant4 10.0 vs. earlier
// 20150112  Remove RM->Initialize() call to allow macro configuration
// 20160111  Remove Geant4 version check since we now hard depend on 10.2+
// 20170816  Add example-specific configuration manager
// 20220718  Remove obsolete pre-processor macros G4VIS_USE and G4UI_USE

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "G4CMPPhysicsList.hh"
#include "G4CMPConfigManager.hh"
#include "Caustic_PhononActionInitialization.hh"
#include "Caustic_PhononConfigManager.hh"
#include "Caustic_PhononDetectorConstruction.hh"

int main(int argc,char** argv)
{
 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

 // Set mandatory initialization classes
 //
 Caustic_PhononDetectorConstruction* detector = new Caustic_PhononDetectorConstruction();
 runManager->SetUserInitialization(detector);

 G4VUserPhysicsList* physics = new G4CMPPhysicsList();
 physics->SetCuts();
 runManager->SetUserInitialization(physics);

 // Set user action classes (different for Geant4 10.0)
 //
 runManager->SetUserInitialization(new Caustic_PhononActionInitialization);

 // Create configuration managers to ensure macro commands exist
 G4CMPConfigManager::Instance();
 Caustic_PhononConfigManager::Instance();

 // Visualization manager
 //
 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();

 // Get the pointer to the User Interface manager
 //--------------Function to Visualisze the Detector Geometry

///////////////////////////////////////////////
// G4UImanager* UImanager = G4UImanager::GetUIpointer();
//
//      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
//   //  UImanager->ApplyCommand("/control/execute vis.mac");
// //  UImanager->ApplyCommand("/run/initialize");
// // UImanager->ApplyCommand("/gps/number 1");
// // UImanager->ApplyCommand("/g4cmp/phononBounces 3000");
// // UImanager->ApplyCommand("/process/inactivate phononScattering");
// // UImanager->ApplyCommand("/process/setVerbose 0 G4CMPSecondaryProduction");
// // UImanager->ApplyCommand("/process/inactivate phononDownconversion");
// // UImanager->ApplyCommand("/gps/energy 0.005 eV");
// // G4UImanager* UImanager = G4UImanager::GetUIpointer();
//
//
//  UImanager->ApplyCommand("/run/initialize");
// UImanager->ApplyCommand("/gps/number 1");
// UImanager->ApplyCommand("/g4cmp/phononBounces 3000");
// UImanager->ApplyCommand("/process/inactivate phononScattering");
// UImanager->ApplyCommand("/process/setVerbose 0 G4CMPSecondaryProduction");
// UImanager->ApplyCommand("/process/inactivate phononDownconversion");
// UImanager->ApplyCommand("/gps/energy 0.005 eV");
//
//
//  UImanager->ApplyCommand("/gps/ang/type iso");
// //  UImanager->ApplyCommand("/gps/direction 0 0 1");
//
//  UImanager->ApplyCommand("/gps/pos/centre 0.0 0.0 0.0 cm");
// // UImanager->ApplyCommand("/gps/pos/centre 0.0 0.0 0.0 cm");
// //
// //
// //
// //
// UImanager->ApplyCommand("/vis/open OGL 600x600-0+0");
// // // UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
// UImanager->ApplyCommand("/vis/drawVolume");
// UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
// UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
// UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTS Red");
// UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTF Green");
// //UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononL Blue");
//  //UImanager->ApplyCommand("/tracking/verbose 2");
//  UImanager->ApplyCommand("/tracking/storeTrajectory 1");
//  UImanager->ApplyCommand("/vis/scene/add/trajectories");
// UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate 4000");
//  UImanager->ApplyCommand("/vis/scene/add/hits");
//
// //UImanager->ApplyCommand("/tracking/verbose 2 ");
// // // UImanager->ApplyCommand("/SiliconSixQubitArrayDetectorConstruction/filmAbsorption 0.5 ");
// // //
// // //
// // // //
// UImanager->ApplyCommand("/vis/scene/add/axes 0 0 0 0.003 m");
// //
// UImanager->ApplyCommand("/run/beamOn 4000");
// // //
// //
// // //Make notes for the meetin g
// //    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
// //    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
//
//
//     ui->SessionStart();




//////////This Part is to run the Macros

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
///////



}
