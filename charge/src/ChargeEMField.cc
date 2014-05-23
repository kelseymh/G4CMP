// $Id$
//
// 20140324  Field is local to specified volume
// 20140331  Do not pass lattice to filed manager; done track-by-track
// 20140522  Migrate to general purpose (in name) field mesh

#include "ChargeEMField.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4CMPFieldManager.hh"
#include "G4LogicalVolume.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  G4ElectricField* fEMfield  = new G4CMPMeshElectricField("Epot_iZip4_small");
  G4FieldManager*  fFieldMgr = new G4CMPFieldManager(fEMfield);

  fFieldMgr->SetDetectorField(fEMfield);
  logVol->SetFieldManager(fFieldMgr, true);
}

ChargeEMField::~ChargeEMField() {;}
