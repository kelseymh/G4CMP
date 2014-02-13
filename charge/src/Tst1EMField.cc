#include "Tst1EMField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "CDMS_iZip4_Field.hh"

#include "G4ios.hh"
#include "XtalFieldManager.hh"


Tst1EMField::Tst1EMField(G4LogicalVolume* logVol){
  CDMS_iZip4_Field*        fEMfield;
  EqEMFieldXtal*          fEquation;
  G4MagIntegratorStepper* fStepper;
  G4FieldManager*         fFieldMgr;
  G4double                fMinStep ;
  G4double                missDistn;
  G4ChordFinder*          fChordFinder ;

  fEMfield = new CDMS_iZip4_Field(G4String("Epot_iZip4"));

  // Create an equation of motion for this field
  fEquation = new EqEMFieldXtal(fEMfield);
  
  G4int nvar = 8;
  fStepper = new G4ClassicalRK4( fEquation, nvar );       

  // Get the global field manager 
  //fFieldMgr= G4TransportationManager::GetTransportationManager()->
  //      GetFieldManager();

  //define accuracy parameters
  //fFieldMgr->SetMinimumEpsilonStep(1e-10*mm);
  //fFieldMgr->SetMaximumEpsilonStep(1e-8*mm);
  //fFieldMgr->SetDeltaOneStep(1e-8*mm);
  //fFieldMgr->SetDeltaIntersection(1e-10*mm);
  fMinStep     = 1e-8*mm ; // minimal step
  //missDistn    = 0.2e-6*mm ; // miss distance

  G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep, 
                                     fStepper, 
                                     fStepper->GetNumberOfVariables() );

  fChordFinder = new G4ChordFinder(fIntgrDriver);

  fFieldMgr = new XtalFieldManager(fEMfield,fChordFinder,true);

  // MHK -- Why is this not being done?
  //logVol->SetFieldManager(fFieldMgr,true);

  G4TransportationManager::GetTransportationManager()->SetFieldManager(fFieldMgr
);
  fFieldMgr->SetDetectorField(fEMfield);
}

Tst1EMField::~Tst1EMField() {;}
