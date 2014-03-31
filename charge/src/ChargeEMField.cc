// $Id$
//
// 20140324  Field is local to specified volume
// 20140331  Do not pass lattice to filed manager; done track-by-track

#include "ChargeEMField.hh"
#include "CDMS_iZip4_Field.hh"
#include "G4CMPFieldManager.hh"
#include "G4LogicalVolume.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  CDMS_iZip4_Field* fEMfield  = new CDMS_iZip4_Field("Epot_iZip4_small");
  G4FieldManager*   fFieldMgr = new G4CMPFieldManager(fEMfield);

  fFieldMgr->SetDetectorField(fEMfield);
  logVol->SetFieldManager(fFieldMgr, true);
}

ChargeEMField::~ChargeEMField() {;}
