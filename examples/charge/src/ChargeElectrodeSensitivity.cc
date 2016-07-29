/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "ChargeElectrodeSensitivity.hh"
#include "ChargeFETDigitizerModule.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPConfigManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

ChargeElectrodeSensitivity::ChargeElectrodeSensitivity(G4String name) :
  G4CMPElectrodeSensitivity(name), fileName(""),
  FET(new ChargeFETDigitizerModule("FETSim"))
{
  SetOutputFile(G4CMPConfigManager::GetHitOutput());
}

ChargeElectrodeSensitivity::~ChargeElectrodeSensitivity()
{
  if (output.is_open()) output.close();
  if (!output.good()) {
    G4ExceptionDescription msg;
    msg << "Error closing output file, " << fileName << ".\n"
        << "Expect bad things like loss of data.";
    G4Exception("ChargeElectrodeSensitivity::~ChargeElectrodSensitivity",
                "Charge004", JustWarning, msg);
  }
  delete FET;
}

void ChargeElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4CMPElectrodeHitsCollection* hitCol =
        static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(GetHCID()));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();
  std::vector<G4CMPElectrodeHit*>::iterator itr = hitVec->begin();
  G4RunManager* runMan = G4RunManager::GetRunManager();
  if (output.good()) {
    for (; itr != hitVec->end(); itr++) {
      output << runMan->GetCurrentRun()->GetRunID() << ','
             << runMan->GetCurrentEvent()->GetEventID() << ','
             << (*itr)->GetTrackID() << ','
             << (*itr)->GetParticleName() << ','
             << (*itr)->GetStartEnergy()/eV << ','
             << (*itr)->GetFinalTime()/ns << ','
             << (*itr)->GetEnergyDeposit()/eV << ','
             << (*itr)->GetStartPosition().getX()/m << ','
             << (*itr)->GetStartPosition().getY()/m << ','
             << (*itr)->GetStartPosition().getZ()/m << ','
             << (*itr)->GetFinalPosition().getX()/m << ','
             << (*itr)->GetFinalPosition().getY()/m << ','
             << (*itr)->GetFinalPosition().getZ()/m << G4endl;
    }
  }
  FET->Digitize();
}

void ChargeElectrodeSensitivity::SetOutputFile(const G4String &fn)
{
  if (fileName != fn) {
    if (output.is_open()) output.close();
    fileName = fn;
    output.open(fileName, std::ios_base::app);
    if (!output.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening output file, " << fileName << ".\n"
          << "Will continue simulation.";
      G4Exception("ChargeElectrodeSensitivity::SetOutputFile", "Charge003",
                  JustWarning, msg);
      output.close();
    } else {
      output << "Run ID,Event ID,Track ID,Particle Name,Start Energy [eV],"
             << "Track Lifetime [ns],Energy Deposit [eV],Start X [m],"
             << "Start Y [m],Start Z [m],End X [m],End Y [m],End Z [m]"
             << G4endl;
    }
  }
}
