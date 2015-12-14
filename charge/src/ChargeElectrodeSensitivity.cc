#include "ChargeElectrodeSensitivity.hh"
#include "ChargeFETDigitizerModule.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPConfigManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

#include <fstream>

ChargeElectrodeSensitivity::ChargeElectrodeSensitivity(G4String name)
:G4CMPElectrodeSensitivity(name), FET(new ChargeFETDigitizerModule("FETSim"))
{}

void ChargeElectrodeSensitivity::Initialize(G4HCofThisEvent *HCE)
{
  //Call base class initialization.
  G4CMPElectrodeSensitivity::Initialize(HCE);

  //Prepare output file.
  output.open(G4CMPConfigManager::GetHitOutput(), std::ios_base::app);
  output << "Run ID,Event ID,Track ID,Charge,Start Energy [eV],Track Lifetime [ns],"
         << "Energy Deposit [eV],Start X [m],Start Y [m],Start Z [m],"
         << "End X [m],End Y [m],End Z [m]"
         << G4endl;

  //Initialize FETSim.
  if (FET->FETSimIsEnabled())
    FET->Initialize();
}

ChargeElectrodeSensitivity::~ChargeElectrodeSensitivity()
{
  output.close();
  delete FET;
}

void ChargeElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4CMPElectrodeHitsCollection* hitCol =
        static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(GetHCID()));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();
  std::vector<G4CMPElectrodeHit*>::iterator itr = hitVec->begin();
  for (; itr != hitVec->end(); itr++)
    output << G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID() << ','
           << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ','
           << (*itr)->GetTrackID() << ',' << (*itr)->GetCharge() << ','
           << (*itr)->GetStartEnergy()/eV << ',' << (*itr)->GetFinalTime()/ns
           << ',' << (*itr)->GetEnergyDeposit()/eV << ','
           << (*itr)->GetStartPosition().getX()/m << ','
           << (*itr)->GetStartPosition().getY()/m << ','
           << (*itr)->GetStartPosition().getZ()/m << ','
           << (*itr)->GetFinalPosition().getX()/m << ','
           << (*itr)->GetFinalPosition().getY()/m << ','
           << (*itr)->GetFinalPosition().getZ()/m << G4endl;
  FET->Digitize();
}

