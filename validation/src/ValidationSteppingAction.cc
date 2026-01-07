/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file ValidationSteppingAction.cc
/// \brief Implementation of the stepping action for the validation example

#include "ValidationSteppingAction.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "ValidationConfigManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include <iostream>
#include <iomanip>

//Default constructor
ValidationSteppingAction::ValidationSteppingAction() {
}

ValidationSteppingAction::~ValidationSteppingAction() {
  fOutputFile.close();
}

//Alternative constructor
void ValidationSteppingAction::UserSteppingAction( const G4Step * step ) {
  //Initialize/open a file here if it has not yet been defined. If I try doing
  //this in the constructor it doesn't seem to inherit the correct name from
  //the config manager...
  if (!fOutputFile.is_open()) {
    G4String stepFileName = ValidationConfigManager::GetStepFileName();
    G4cout << "StepFileName: " << stepFileName << G4endl;
    
    //Upon construction of this class, create a txt file with step information 
    fOutputFile.open(stepFileName.c_str(),std::ios::trunc);
  }
  
  //First up: do generic exporting of step information (no cuts made here)
  ExportStepInformation(step);

  clock_t timestamp;
  timestamp = clock();
  
  return;
}

// Do a set of queries of information to test for anharmonic decay
void ValidationSteppingAction::ExportStepInformation( const G4Step * step ) {
  //Test
  G4StepPoint * preSP = step->GetPreStepPoint();
  G4StepPoint * postSP = step->GetPostStepPoint();

  int runNo = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  int eventNo = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  
  if (eventNo > 1000000) { return; }
  int trackNo = step->GetTrack()->GetTrackID();
  std::string particleName 
    = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  double preStepX_mm = preSP->GetPosition().x() / CLHEP::mm;
  double preStepY_mm = preSP->GetPosition().y() / CLHEP::mm;
  double preStepZ_mm = preSP->GetPosition().z() / CLHEP::mm;
  double preStepT_ns = preSP->GetGlobalTime() / CLHEP::ns;
  double preStepEnergy_eV = preSP->GetTotalEnergy() / CLHEP::eV;
  double preStepKinEnergy_eV = preSP->GetKineticEnergy() / CLHEP::eV;

  
  double postStepX_mm = postSP->GetPosition().x() / CLHEP::mm;
  double postStepY_mm = postSP->GetPosition().y() / CLHEP::mm;
  double postStepZ_mm = postSP->GetPosition().z() / CLHEP::mm;
  double postStepT_ns = postSP->GetGlobalTime() / CLHEP::ns;
  double postStepEnergy_eV = postSP->GetTotalEnergy() / CLHEP::eV;
  double postStepKinEnergy_eV = postSP->GetKineticEnergy() / CLHEP::eV;

  //Get reflection count
  size_t nReflections
    = G4CMP::GetTrackInfo<G4CMPVTrackInfo>(step->GetTrack())->ReflectionCount();
    
  std::string stepProcess = postSP->GetProcessDefinedStep()->GetProcessName();

  //Fill the output file with the step info  
  fOutputFile << runNo << " " << eventNo << " " << trackNo << " "
              << particleName << " " << std::setprecision(14) << preStepX_mm
              << " " << preStepY_mm << " " << preStepZ_mm << " "
              << preStepT_ns << " " << preStepEnergy_eV << " "
              << preStepKinEnergy_eV << " " << postStepX_mm << " "
              << postStepY_mm << " " << postStepZ_mm << " " << postStepT_ns
              << " " << postStepEnergy_eV << " " << postStepKinEnergy_eV
              << " " << nReflections << " " << stepProcess << std::endl;
}
