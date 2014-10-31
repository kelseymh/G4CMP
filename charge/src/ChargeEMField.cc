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
#include "G4UniformElectricField.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "ChargeFieldMessenger.hh"

ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) :
  EPotFile(G4CMPConfigManager::GetEpotFile()),
  UseConstEField(G4CMPConfigManager::GetVoltage()!=0.),
  ConstEFieldMag(G4CMPConfigManager::GetVoltage()/(25.4*mm)),
  ConstEFieldDir(0,0,1), DetectorVolume(logVol) {
  Messenger = new ChargeFieldMessenger(this);
  Build();			// Ensure that default field is alway created
}

ChargeEMField::~ChargeEMField() {;}

void ChargeEMField::Build() {
  G4ElectricField* fEMfield = 0;

  if (UseConstEField) {
    fEMfield = new G4UniformElectricField(ConstEFieldMag*ConstEFieldDir);
  } else {
    fEMfield = new G4CMPMeshElectricField(EPotFile);
  }

  // Ensure that logical volume has a field manager attached
  if (!DetectorVolume->GetFieldManager()) {
    G4FieldManager* fFieldMgr = new G4CMPFieldManager(fEMfield);
    DetectorVolume->SetFieldManager(fFieldMgr, true);
  }

  DetectorVolume->GetFieldManager()->SetDetectorField(fEMfield);
}
