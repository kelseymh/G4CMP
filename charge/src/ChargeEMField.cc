// $Id$

#include "ChargeEMField.hh"
#include "CDMS_iZip4_Field.hh"
#include "G4CMPFieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  G4LatticeLogical* lattice = G4LatticeManager::GetLatticeManager()->GetLattice(logVol->GetMaterial());

  CDMS_iZip4_Field* fEMfield  = new CDMS_iZip4_Field(G4String("Epot_iZip4"));
  G4FieldManager*   fFieldMgr = new G4CMPFieldManager(fEMfield, lattice);

  fFieldMgr->SetDetectorField(fEMfield);

  // MHK -- Why is this not being done?
  //logVol->SetFieldManager(fFieldMgr,true);

  G4TransportationManager::GetTransportationManager()->SetFieldManager(fFieldMgr
);
}

ChargeEMField::~ChargeEMField() {;}
