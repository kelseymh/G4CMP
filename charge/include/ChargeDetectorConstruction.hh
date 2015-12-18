// $Id: 165a1df2a9364e285f19b48f5c016259faee053b $

#ifndef ChargeDetectorConstruction_h
#define ChargeDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class ChargeElectrodeSensitivity;
class G4LatticeManager;
class G4Material;
class G4VPhysicalVolume;
class G4ElectricField;


class ChargeDetectorConstruction : public G4VUserDetectorConstruction {
public:
  ChargeDetectorConstruction();
  virtual ~ChargeDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
     
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachField(G4LogicalVolume* lv);
  void AttachLattice(G4VPhysicalVolume* pv);
  void AttachSensitivity(G4LogicalVolume* lv);

private:
  G4LatticeManager* latManager;
  G4ElectricField* fEMField;
  ChargeElectrodeSensitivity* sensitivity;
  G4Material* liquidHelium;
  G4Material* germanium;
  G4Material* aluminum;
  G4Material* tungsten;
  G4VPhysicalVolume* worldPhys;
  G4double zipThickness; // Useful for geom. and field
  G4double epotScale;
  G4double voltage;
  G4bool constructed;
  G4String epotFileName;
  G4String outputFileName;
};

#endif

