/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

//20240110 Israel Hernandez -- Illinois Institute of Technology, Quantum Science Center and Fermilab
#ifndef Caustic_PhononDetectorConstruction_h
#define Caustic_PhononDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;


class Caustic_PhononDetectorConstruction : public G4VUserDetectorConstruction {
public:
  Caustic_PhononDetectorConstruction();
  virtual ~Caustic_PhononDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();

private:
  void Caustic_DefineMaterials();
  void Caustic_SetupGeometry();


private:
  G4Material* fLiquidHelium;
  G4Material* fBolometer;
  G4Material* fAluminum;
  G4Material* fOxigen;
  G4Material* fSubstrate;
  G4VPhysicalVolume* fWorldPhys;
  G4CMPSurfaceProperty* topSurfProp;
  G4CMPSurfaceProperty* wallSurfProp;
  G4CMPElectrodeSensitivity* electrodeSensitivity;

  G4bool fConstructed;		// Flag to not re-recreate surface properties
};

#endif
