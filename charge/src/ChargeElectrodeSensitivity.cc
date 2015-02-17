#include "ChargeElectrodeSensitivity.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPConfigManager.hh"
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
:G4CMPElectrodeSensitivity(name)
{
  output.open(G4CMPConfigManager::GetHitOutput());
  output << "Track ID, Charge, Start Energy [eV], Track Lifetime [ns], "
         << "Energy Deposit [eV], Start Position [m], End Position [m]"
         << G4endl;
}

ChargeElectrodeSensitivity::~ChargeElectrodeSensitivity()
{
  output.close();
}

void ChargeElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4CMPElectrodeHitsCollection* hitCol =
        static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(GetHCID()));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();
  /* Can we use C++11 yet? Look how nice this should look:
  for (auto j : *hitVec)
    output << j->GetTrackID() << ", " << j->GetCharge() << ", "
           << j->GetStartEnergy()/eV << ", " << j->GetTrackTime()/ns
           << ", " << j->GetEnergyDeposit()/eV << ", "
           << j->GetStartPosition().getX()/m << ", "
           << j->GetStartPosition().getY()/m << ", "
           << j->GetStartPosition().getZ()/m << ", "
           << j->GetFinalPosition().getX()/m << ", "
           << j->GetFinalPosition().getY()/m << ", "
           << j->GetFinalPosition().getZ()/m << G4endl;
  */
  std::vector<G4CMPElectrodeHit*>::iterator itr = hitVec->begin();
  for (; itr != hitVec->end(); itr++)
    output << (*itr)->GetTrackID() << ", " << (*itr)->GetCharge() << ", "
           << (*itr)->GetStartEnergy()/eV << ", " << (*itr)->GetTrackTime()/ns
           << ", " << (*itr)->GetEnergyDeposit()/eV << ", "
           << (*itr)->GetStartPosition().getX()/m << ", "
           << (*itr)->GetStartPosition().getY()/m << ", "
           << (*itr)->GetStartPosition().getZ()/m << ", "
           << (*itr)->GetFinalPosition().getX()/m << ", "
           << (*itr)->GetFinalPosition().getY()/m << ", "
           << (*itr)->GetFinalPosition().getZ()/m << G4endl;
}

