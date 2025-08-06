/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20241024 Israel Hernandez -- IIT, QSC and Fermilab
// 20250101 Michael Kelsey -- Instantiate SD in ConstructSDandField()

#ifndef Caustic_PhononDetectorConstruction_h
#define Caustic_PhononDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;


class Caustic_PhononDetectorConstruction : public G4VUserDetectorConstruction {
public:
  Caustic_PhononDetectorConstruction();
  virtual ~Caustic_PhononDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

private:
  void Caustic_DefineMaterials();
  void Caustic_SetupGeometry();

private:
  G4Material* fLiquidHelium;
  G4Material* fBolometer;
  G4Material* fCrystalMaterial;
  G4VPhysicalVolume* fWorldPhys;
  G4LogicalVolume* fpSubstrateLV;
  G4CMPSurfaceProperty* topSurfProp;
  G4CMPSurfaceProperty* wallSurfProp;

  G4bool fConstructed;		// Flag to not re-recreate surface properties
};

#endif
