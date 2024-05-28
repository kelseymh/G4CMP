/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


// Basic User Stepping action for the silicon six qubit array (mostly for debugging)

#include "RISQTutorialSteppingAction.hh"
#include <iostream>
#include "globals.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Threading.hh"

#include "G4RunManager.hh"
#include "G4StepPoint.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Default constructor
RISQTutorialSteppingAction::RISQTutorialSteppingAction()
{
  //Upon construction of this class, create a ROOT file with step information and a tree with variables for
  //storing the step information if needed
  //fOutputFile.open("StepInformationFile.txt",std::ios::trunc);

  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
RISQTutorialSteppingAction::~RISQTutorialSteppingAction()
{
  //fOutputFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Alternative constructor
void RISQTutorialSteppingAction::UserSteppingAction( const G4Step * step )
{
  //For now, simple: look at the pre-step point volume name and the track name
  //  std::cout << "REL stepping. PreSP volume name: " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << ", track particle type: " << step->GetTrack()->GetParticleDefinition()->GetParticleName() << std::endl;

  //First up: do generic exporting of step information (no cuts made here)
  //ExportStepInformation(step);
  
  
  return;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Do a set of queries of information to test for anharmonic decay
void RISQTutorialSteppingAction::ExportStepInformation( const G4Step * step )
{
  //Test
  G4StepPoint * preSP = step->GetPreStepPoint();
  G4StepPoint * postSP = step->GetPostStepPoint();

  int runNo = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  int eventNo = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  int trackNo = step->GetTrack()->GetTrackID();
  std::string particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  double preStepX_mm = preSP->GetPosition().x() / CLHEP::mm;
  double preStepY_mm = preSP->GetPosition().y() / CLHEP::mm;
  double preStepZ_mm = preSP->GetPosition().z() / CLHEP::mm;
  double preStepEnergy_eV = preSP->GetTotalEnergy() / CLHEP::eV;
  double preStepKinEnergy_eV = preSP->GetKineticEnergy() / CLHEP::eV;
  
  double postStepX_mm = postSP->GetPosition().x() / CLHEP::mm;
  double postStepY_mm = postSP->GetPosition().y() / CLHEP::mm;
  double postStepZ_mm = postSP->GetPosition().z() / CLHEP::mm; 
  double postStepEnergy_eV = postSP->GetTotalEnergy() / CLHEP::eV;
  double postStepKinEnergy_eV = postSP->GetKineticEnergy() / CLHEP::eV;

  std::string stepProcess = postSP->GetProcessDefinedStep()->GetProcessName();
  
  //Fill the output file with the step info
  fOutputFile << runNo << " " << eventNo << " " << trackNo << " " << particleName << " "
	      << preStepX_mm << " " << preStepY_mm << " " << preStepZ_mm << " " << preStepEnergy_eV << " " << preStepKinEnergy_eV << " " 
	      << postStepX_mm << " " << postStepY_mm << " " << postStepZ_mm << " " << postStepEnergy_eV << " " << postStepKinEnergy_eV
	      << " " << stepProcess << std::endl;

  
  
}
