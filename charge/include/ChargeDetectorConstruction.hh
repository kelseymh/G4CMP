
#ifndef ChargeDetectorConstruction_h
#define ChargeDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LatticeManager.hh"
#include "ChargeEMField.hh"

class G4Material;
class G4VPhysicalVolume;
class ChargeDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ChargeDetectorConstruction();
    virtual ~ChargeDetectorConstruction();

public:
    virtual G4VPhysicalVolume* Construct();
     
private:
  void DefineMaterials();
  void SetupGeometry();

  private:


    G4Material* liquidHelium;
    G4Material* germanium;
    G4Material* alminum;
    G4Material* tungsten;
    G4VPhysicalVolume* worldPhys;
    G4bool constructed;
    G4bool ifField;
    ChargeEMField* field;

  public:
    inline void Field(G4bool bl) { ifField = bl; }
};

#endif

