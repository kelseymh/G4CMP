// $Id$
//
// 20140324  Field is local to specified volume
// 20140331  Do not pass lattice to filed manager; done track-by-track
// 20140522  Migrate to general purpose (in name) field mesh
// 20140623  Add local variable to set name string from envvar

#include "ChargeEMField.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4CMPFieldManager.hh"
#include "G4LogicalVolume.hh"

// Get name of field file from envvar
namespace {
  const G4String epotFile =
    (getenv("G4CMP_EPOT_FILE") ? getenv("G4CMP_EPOT_FILE") : "Epot_iZip4_small");
}

ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  G4ElectricField* fEMfield  = new G4CMPMeshElectricField(epotFile);
  G4FieldManager*  fFieldMgr = new G4CMPFieldManager(fEMfield);

  fFieldMgr->SetDetectorField(fEMfield);
  logVol->SetFieldManager(fFieldMgr, true);
}

ChargeEMField::~ChargeEMField() {;}
