/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


// Basic User Stepping action for Quasiparticle

#include "QuasiparticleSteppingAction.hh"
#include <iostream>
#include <iomanip>
#include "globals.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Threading.hh"

#include "G4RunManager.hh"
#include "G4StepPoint.hh"
#include "G4CMPVTrackInfo.hh"
#include "G4CMPTrackUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Default constructor
QuasiparticleSteppingAction::QuasiparticleSteppingAction()
{
  //Upon construction of this class, create a txt file with step information 
  fOutputFile.open("StepInformationFile.txt",std::ios::trunc);

  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
QuasiparticleSteppingAction::~QuasiparticleSteppingAction()
{
  fOutputFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Alternative constructor
void QuasiparticleSteppingAction::UserSteppingAction( const G4Step * step )
{
  //For now, simple: look at the pre-step point volume name and the track name
  //  std::cout << "REL stepping. PreSP volume name: " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << ", track particle type: " << step->GetTrack()->GetParticleDefinition()->GetParticleName() << std::endl;

  //First up: do generic exporting of step information (no cuts made here)
  ExportStepInformation(step);

  clock_t timestamp;
  timestamp = clock();
  //G4cout << "-----------------------> Time in user stepping action: " << double(timestamp)/double(CLOCKS_PER_SEC) << " seconds" << G4endl;
  
  return;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do a set of queries of information to test for anharmonic decay
void QuasiparticleSteppingAction::ExportStepInformation( const G4Step * step )
{
  //Test
  G4StepPoint * preSP = step->GetPreStepPoint();
  G4StepPoint * postSP = step->GetPostStepPoint();

  int runNo = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  int eventNo = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  
  if( eventNo > 1000000 ){ return; }
  int trackNo = step->GetTrack()->GetTrackID();
  std::string particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  double preStepX_mm = preSP->GetPosition().x() / CLHEP::mm;
  double preStepY_mm = preSP->GetPosition().y() / CLHEP::mm;
  double preStepZ_mm = preSP->GetPosition().z() / CLHEP::mm;
  double preStepT_ns = preSP->GetGlobalTime() / CLHEP::ns;
  double preStepEnergy_eV = preSP->GetTotalEnergy() / CLHEP::eV;
  double preStepKinEnergy_eV = preSP->GetKineticEnergy() / CLHEP::eV;

  //std::cout << "Pre-step energy: " << preStepKinEnergy_eV << " eV" << std::endl;
  
  double postStepX_mm = postSP->GetPosition().x() / CLHEP::mm;
  double postStepY_mm = postSP->GetPosition().y() / CLHEP::mm;
  double postStepZ_mm = postSP->GetPosition().z() / CLHEP::mm;
  double postStepT_ns = postSP->GetGlobalTime() / CLHEP::ns;
  double postStepEnergy_eV = postSP->GetTotalEnergy() / CLHEP::eV;
  double postStepKinEnergy_eV = postSP->GetKineticEnergy() / CLHEP::eV;

  //Get reflection count
  size_t nReflections = G4CMP::GetTrackInfo<G4CMPVTrackInfo>(step->GetTrack())->ReflectionCount();
    
  std::string stepProcess = postSP->GetProcessDefinedStep()->GetProcessName();


  //Fill the output file with the step info  
  fOutputFile << runNo << " " << eventNo << " " << trackNo << " " << particleName << " "
	      << std::setprecision(14) << preStepX_mm << " " << preStepY_mm << " " << preStepZ_mm << " " << preStepT_ns << " " << preStepEnergy_eV << " "
	      << preStepKinEnergy_eV << " " << postStepX_mm << " " << postStepY_mm << " " << postStepZ_mm << " " << postStepT_ns << " "
	      << postStepEnergy_eV << " " << postStepKinEnergy_eV << " " << nReflections << " " << stepProcess << std::endl;
  
  
  
}
