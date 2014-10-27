// $Id$
//
// 20140324  Field is local to specified volume
// 20140331  Do not pass lattice to filed manager; done track-by-track
// 20140522  Migrate to general purpose (in name) field mesh
// 20140623  Add local variable to set name string from envvar

#include "ChargeEMField.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4CMPFieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "ChargeFieldMessenger.hh"

ChargeEMField::ChargeEMField(G4LogicalVolume* logVol) :
                              UseConstEField(false),
                              ConstEFieldMag(158*volt/m),
                              ConstEFieldDir(G4ThreeVector(0,0,1)),
                              DetectorVolume(logVol) {
  EPotFile =
    (getenv("G4CMP_EPOT_FILE") ? getenv("G4CMP_EPOT_FILE") : "Epot_iZip4");
  Messenger = new ChargeFieldMessenger(this);
}

ChargeEMField::~ChargeEMField() {;}

void ChargeEMField::Build() {
  G4ElectricField* fEMfield;

  if (UseConstEField) {
    fEMfield = new G4UniformElectricField(ConstEFieldMag*ConstEFieldDir);
  } else {
    fEMfield = new G4CMPMeshElectricField(EPotFile);
  }

  if(DetectorVolume->GetFieldManager()) {
    DetectorVolume->GetFieldManager()->SetDetectorField(fEMfield);
  } else {
    G4FieldManager*  fFieldMgr = new G4CMPFieldManager(fEMfield);
    fFieldMgr->SetDetectorField(fEMfield);
    DetectorVolume->SetFieldManager(fFieldMgr, true);
  }
}
