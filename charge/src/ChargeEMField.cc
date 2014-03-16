#include "ChargeEMField.hh"
#include "CDMS_iZip4_Field.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4CMPFieldManager.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  CDMS_iZip4_Field* fEMfield  = new CDMS_iZip4_Field(G4String("Epot_iZip4"));
  G4FieldManager*   fFieldMgr = new G4CMPFieldManager(fEMfield);

  // MHK -- Why is this not being done?
  //logVol->SetFieldManager(fFieldMgr,true);

  G4TransportationManager::GetTransportationManager()->SetFieldManager(fFieldMgr
);
  fFieldMgr->SetDetectorField(fEMfield);
}

ChargeEMField::~ChargeEMField() {;}
