// $Id$
//
// 20140324  Field is local to specified volume
// 20140331  Do not pass lattice to filed manager; done track-by-track
// 20140522  Migrate to general purpose (in name) field mesh
// 20140623  Add local variable to set name string from envvar
// 20141029  Get file name or fixed voltage from configuration manager

#include "ChargeEMField.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4CMPFieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UniformElectricField.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  G4double voltage = G4CMPConfigManager::GetVoltage();
  const G4String& epotFile = G4CMPConfigManager::GetEpotFile();

  G4ElectricField* fEMfield = 0;
  if (voltage == 0.) fEMfield = new G4CMPMeshElectricField(epotFile);
  else {
    G4ThreeVector evec(0.,0.,voltage/(25.4*mm));	// CDMS iZIP thickness
    fEMfield = new G4UniformElectricField(evec);
  }

  G4FieldManager* fFieldMgr = new G4CMPFieldManager(fEMfield);

  fFieldMgr->SetDetectorField(fEMfield);
  logVol->SetFieldManager(fFieldMgr, true);
}

ChargeEMField::~ChargeEMField() {;}
