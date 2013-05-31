
#ifndef Tst1DetectorConstruction_h
#define Tst1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "LatticeManager2.hh"
#include "Tst1EMField.hh"

class G4Material;
class G4VPhysicalVolume;
class Tst1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst1DetectorConstruction();
    virtual ~Tst1DetectorConstruction();

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
    Tst1EMField* field;

  public:
    inline void Field(G4bool bl) { ifField = bl; }
};

#endif

