/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructureDetectorConstruction.hh
/// \brief Definition of the HeterostructureDetectorConstruction class

// 20260815  Selby Q. Dang (Stanford/SLAC)

#ifndef HeterostructureDetectorConstruction_h
#define HeterostructureDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;


class HeterostructureDetectorConstruction : public G4VUserDetectorConstruction {
public:
  HeterostructureDetectorConstruction();
  virtual ~HeterostructureDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
private:
  void DefineMaterials();
  void SetupGeometry();

private:
  G4Material* fLiquidHelium;
  G4Material* fGermanium;
  G4Material* fSilicon;
  G4VPhysicalVolume* fWorldPhys;
  G4CMPSurfaceProperty* fGeSiInterface;
  G4CMPSurfaceProperty* fSiGeInterface;
  G4CMPSurfaceProperty* fGeVacInterface;
  G4CMPSurfaceProperty* fSiVacInterface;
  G4bool fConstructed;		// Flag to not re-recreate surface properties
};

#endif

