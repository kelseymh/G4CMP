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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"
//#include "B1Analysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"

#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1100
#include "G4AnalysisManager.hh"
#else // G4VERSION_NUMBER >= 1100
#include "g4root.hh"
// #include "g4csv.hh"
#endif // G4VERSION_NUMBER >= 1100

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int B1EventAction::fShotNumber = 0; // Initialize static variable

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fPhononsInterfaceWeighted(0.),
  fPhononEnergyInterfaceWeighted(0.),
  fPhononsAbove2DeltaWeighted(0.),
  fPhononEnergyAbove2DeltaWeighted(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  fShotNumber++; // Increment the shot number
  
  //Start G4CMP application
  // Clear per-event electron bookkeeping
  fElectronTracksInSi.clear();
  fElectronTracksReachedInterface.clear();
  fElectronTracksTrappedInSi.clear();

  fElectronInterfaceKE.clear();
  fElectronInterfaceVperp.clear();
  fElectronInterfaceTime.clear();
  fElectronInterfacePos.clear();

  fPhononsInterfaceWeighted = 0.;
  fPhononEnergyInterfaceWeighted = 0.;
  fPhononsAbove2DeltaWeighted = 0.;
  fPhononEnergyAbove2DeltaWeighted = 0.;
  //End G4CMP application
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{
  fRunAction->AddEdep(fEdep);

  G4double nInSi    = GetNElectronsInSi();
  G4double nReached = GetNElectronsReachedInterface();
  G4double nTrapped = GetNElectronsTrappedInSi();
  G4double fracHit  = GetElectronReachFraction();
  G4int eventID     = event->GetEventID();

  auto analysisManager = G4AnalysisManager::Instance();

  // Electron event summary
  analysisManager->FillNtupleIColumn(1, 0, eventID);
  analysisManager->FillNtupleDColumn(1, 1, nInSi);
  analysisManager->FillNtupleDColumn(1, 2, nReached);
  analysisManager->FillNtupleDColumn(1, 3, nTrapped);
  analysisManager->FillNtupleDColumn(1, 4, fracHit);
  analysisManager->AddNtupleRow(1);

  // Phonon event summary
  analysisManager->FillNtupleIColumn(2, 0, eventID);
  analysisManager->FillNtupleDColumn(2, 1, fPhononsInterfaceWeighted);
  analysisManager->FillNtupleDColumn(2, 2, fPhononEnergyInterfaceWeighted / eV);
  analysisManager->FillNtupleDColumn(2, 3, fPhononsAbove2DeltaWeighted);
  analysisManager->FillNtupleDColumn(2, 4, fPhononEnergyAbove2DeltaWeighted / eV);
  analysisManager->AddNtupleRow(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::CountElectronInSi(G4int trackID, G4double weight)
{
  if (fElectronTracksInSi.find(trackID) == fElectronTracksInSi.end()) {
    fElectronTracksInSi[trackID] = weight;
  }
}

void B1EventAction::CountElectronReachedInterface(G4int trackID, G4double weight)
{
  if (fElectronTracksReachedInterface.find(trackID) == fElectronTracksReachedInterface.end()) {
    fElectronTracksReachedInterface[trackID] = weight;
  }

  auto it = fElectronTracksTrappedInSi.find(trackID);
  if (it != fElectronTracksTrappedInSi.end()) {
    fElectronTracksTrappedInSi.erase(it);
  }
}

void B1EventAction::CountElectronTrappedInSi(G4int trackID, G4double weight)
{
  if (fElectronTracksReachedInterface.find(trackID) == fElectronTracksReachedInterface.end()) {
    if (fElectronTracksTrappedInSi.find(trackID) == fElectronTracksTrappedInSi.end()) {
      fElectronTracksTrappedInSi[trackID] = weight;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double B1EventAction::GetNElectronsInSi() const
{
  G4double sum = 0.0;
  for (const auto& kv : fElectronTracksInSi) sum += kv.second;
  return sum;
}

G4double B1EventAction::GetNElectronsReachedInterface() const
{
  G4double sum = 0.0;
  for (const auto& kv : fElectronTracksReachedInterface) sum += kv.second;
  return sum;
}

G4double B1EventAction::GetNElectronsTrappedInSi() const
{
  G4double sum = 0.0;
  for (const auto& kv : fElectronTracksTrappedInSi) sum += kv.second;
  return sum;
}

G4double B1EventAction::GetElectronReachFraction() const
{
  G4double nInSi = GetNElectronsInSi();
  if (nInSi <= 0.0) return 0.0;

  return GetNElectronsReachedInterface() / nInSi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::AddInterfaceElectronKE(G4double ke)
{
  fElectronInterfaceKE.push_back(ke);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::AddInterfaceElectronVperp(G4double vperp)
{
  fElectronInterfaceVperp.push_back(vperp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::AddInterfaceElectronTime(G4double time)
{
  fElectronInterfaceTime.push_back(time);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::AddInterfaceElectronPosition(const G4ThreeVector& pos)
{
  fElectronInterfacePos.push_back(pos);
}
              
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::AddPhononInterfaceHit(G4double weight, G4double energy)
{
  const G4double twoDelta = 3.6e-4 * eV;

  fPhononsInterfaceWeighted += weight;
  fPhononEnergyInterfaceWeighted += weight * energy;

  if (energy >= twoDelta) {
    fPhononsAbove2DeltaWeighted += weight;
    fPhononEnergyAbove2DeltaWeighted += weight * energy;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
