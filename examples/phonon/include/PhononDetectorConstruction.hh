/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file PhononDetectorConstruction.hh
/// \brief Definition of the PhononDetectorConstruction class
//
// 20221006  M. Kelsey -- Remove "IsField" flag, unnecessary with phonons.
//		Add material properties for aluminum phonon sensors
// 20260827  B. Zatschler -- Add new function ConstructSDandField().
//              Remove data member electrodeSensitivity.

#ifndef PhononDetectorConstruction_h
#define PhononDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;


class PhononDetectorConstruction : public G4VUserDetectorConstruction {
public:
  PhononDetectorConstruction();
  virtual ~PhononDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty* surfProp);

private:
  G4Material* fLiquidHelium;
  G4Material* fGermanium;
  G4Material* fAluminum;
  G4Material* fTungsten;
  G4VPhysicalVolume* fWorldPhys;
  G4CMPSurfaceProperty* topSurfProp;
  G4CMPSurfaceProperty* botSurfProp;
  G4CMPSurfaceProperty* wallSurfProp;

  G4bool fConstructed;		// Flag to not re-recreate surface properties
};

#endif  /* PhononDetectorConstruction_h */
