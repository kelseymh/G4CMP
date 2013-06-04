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

#include "XWrapperProcess.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "G4SystemOfUnits.hh"


XWrapperProcess::XWrapperProcess(const G4String& aName)
:G4VDiscreteProcess(aName)
{   
  
  channelingFactor = 20.0;
  
    if (verboseLevel>1) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperProcess::XWrapperProcess(const G4String& processName, G4VDiscreteProcess* toRegister)
  :G4VDiscreteProcess("XWrapperProcess")
{
  registeredProcess = toRegister;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XWrapperProcess::~XWrapperProcess()
{

}

XWrapperProcess::XWrapperProcess(XWrapperProcess& right)
: G4VDiscreteProcess(right)
{;}

void XWrapperProcess::StartTracking(G4Track* aTrack)
{
  registeredProcess->StartTracking(aTrack);
  currentInteractionLength = -1.0;
  theNumberOfInteractionLengthLeft = -1.0;
  theInitialNumberOfInteractionLength=-1.0;  
}

void XWrapperProcess::registerProcess(G4VDiscreteProcess* toRegister){
  registeredProcess = toRegister;
}

G4double
XWrapperProcess::GetMeanFreePath(
			    const G4Track& aTrack, 
			    G4double prevStepSize, 
			    G4ForceCondition* condition  )
{

  return DBL_MAX;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XWrapperProcess::PostStepGetPhysicalInteractionLength (const G4Track &aTrack, 
								G4double previousStepSize, 
								G4ForceCondition *condition)
{

  if ( (previousStepSize < 0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
    // beginning of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else if ( previousStepSize > 0.0) {
    // subtract NumberOfInteractionLengthLeft 
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
  } else {
    // zero step
    //  DO NOTHING
  }

  //determine factor by which to change mfp
  ChannelingParticleUserInfo* chanInfo = (ChannelingParticleUserInfo*) aTrack.GetUserInformation();
  channelingFactor = chanInfo->GetChannelingFactor();
  G4double preStepChannelingFactor = chanInfo->GetPreStepChannelingFactor();
  ///////////////////////////
  //for debug purposes only: set channeling factor by hand
  ///////////////////////////
  channelingFactor = 1.0;  // <--REMOVE THIS AFTER DEBUGGING!!!!!
  ///////////////////////////

  G4double regIntLength = registeredProcess->PostStepGetPhysicalInteractionLength(aTrack, previousStepSize/preStepChannelingFactor, condition);
  G4double regIntNumber = registeredProcess->GetNumberOfInteractionLengthLeft();
  currentInteractionLength = regIntLength / regIntNumber;
  theNumberOfInteractionLengthLeft = regIntNumber;

  G4cout<<"\nXRawpperProcess::PostStepGetPhysicalInteractionLength:\n\tregistered PIL: \t"<<regIntLength;
  G4cout<<"\n\tregistered number of  Int Length Left:\t"<<regIntNumber;
  G4cout<<"\n\tcurrent int length:\t"<<currentInteractionLength/m;
  G4cout<<"\n\tprevious step size:\t"<<previousStepSize/m;
  G4cout<<"\n";

  currentInteractionLength = channelingFactor * theNumberOfInteractionLengthLeft * currentInteractionLength;
  return currentInteractionLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
XWrapperProcess::PostStepDoIt( const G4Track& aTrack,
                                const G4Step& aStep )
{
 
    
  //aParticleChange.Initialize(aTrack);

  //G4VParticleChange* buffer = registeredProcess->PostStepDoIt(aTrack, aStep);

  //registeredProcess->ResetNumberOfInteractionLengthLeft();
  //  ClearNumberOfInteractionLengthLeft();
    
  return registeredProcess->PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool XWrapperProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
    return(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XWrapperProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

