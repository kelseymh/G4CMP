/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQTutorialDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//

#ifndef RISQTutorialDetectorConstruction_h
#define RISQTutorialDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class RISQTutorialDetectorConstruction : public G4VUserDetectorConstruction {
public:
  RISQTutorialDetectorConstruction();
  virtual ~RISQTutorialDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty * surfProp);

  
private:
  G4Material* fLiquidHelium;
  G4Material* fGermanium;
  G4Material* fAluminum;
  G4Material* fTungsten;
  G4Material* fSilicon;
  G4Material* fNiobium;
  G4VPhysicalVolume* fWorldPhys;

  G4CMPSurfaceProperty* fSiNbInterface;
  G4CMPSurfaceProperty* fSiCopperInterface;
  G4CMPSurfaceProperty* fSiVacuumInterface;

  
  G4CMPElectrodeSensitivity* fSuperconductorSensitivity;
  G4bool fConstructed;
  //G4bool fIfField;
  
  //public:
  //inline void Field(G4bool bl) { fIfField = bl; }
};

#endif

