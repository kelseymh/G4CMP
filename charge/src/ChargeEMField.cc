// $Id$
//
// 20130324  Field is local to specified volume

#include "ChargeEMField.hh"
#include "CDMS_iZip4_Field.hh"
#include "G4CMPFieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LogicalVolume.hh"


ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) {
  G4LatticeLogical* lattice = G4LatticeManager::GetLatticeManager()->GetLattice(logVol->GetMaterial());

  CDMS_iZip4_Field* fEMfield  = new CDMS_iZip4_Field("Epot_iZip4_small");
  G4FieldManager*   fFieldMgr = new G4CMPFieldManager(fEMfield, lattice);

  fFieldMgr->SetDetectorField(fEMfield);
  logVol->SetFieldManager(fFieldMgr, true);
}

ChargeEMField::~ChargeEMField() {;}
