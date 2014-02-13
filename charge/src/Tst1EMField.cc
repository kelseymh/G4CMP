#include "Tst1EMField.hh"
#include "G4UniformMagField.hh"
#include "G4UniformElectricField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4LogicalVolume.hh"
#include "G4ExplicitEuler.hh" 
#include "G4CashKarpRKF45.hh"
#include "G4SimpleHeum.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"
#include "XtalFieldManager.hh"


Tst1EMField::Tst1EMField(G4LogicalVolume* logVol){

  G4ElectricField*        fEMfield;
  EqEMFieldXtal*          fEquation;
  //G4EqMagElectricField*   fEquation;
  G4MagIntegratorStepper* fStepper;
  G4FieldManager*         fFieldMgr;
  //XtalFieldManager*       fFieldMgr;
  G4double                fMinStep ;
  G4double                missDistn;
  G4ChordFinder*          fChordFinder ;


  //fEMfield = new G4UniformElectricField(
  //               G4ThreeVector(0.0,0.0, 5.0*volt/cm));

    fEMfield = new G4UniformElectricField(
                 G4ThreeVector(0.0,0.0, 0.524*volt/cm));

  
  // Create an equation of motion for this field
  fEquation = new EqEMFieldXtal(fEMfield);
  //fEquation = new G4EqMagElectricField(fEMfield); 
  
  G4int nvar = 8;
  fStepper = new DMCClassicalRK4( fEquation, nvar );       
  //fStepper = new G4ClassicalRK4( fEquation, nvar );       

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


  // Set this field to the global field manager 
  //fFieldMgr->SetDetectorField(fEMfield );
  

  
  G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep, 
                                     fStepper, 
                                     fStepper->GetNumberOfVariables() );

  fChordFinder = new G4ChordFinder(fIntgrDriver);
  //fChordFinder->SetDeltaChord(missDistn);
  fFieldMgr = new XtalFieldManager(fEMfield,fChordFinder,true);
  //logVol->SetFieldManager(fFieldMgr,true);

  
  
G4TransportationManager::GetTransportationManager()->SetFieldManager(fFieldMgr
) ;
  fFieldMgr->SetDetectorField(fEMfield );
//   G4TransportationManager::GetTransportationManager()->SetFieldManager(new 
//XtalFieldManager(fEMfield, fChordFinder, true));
}


Tst1EMField::~Tst1EMField(){

}
