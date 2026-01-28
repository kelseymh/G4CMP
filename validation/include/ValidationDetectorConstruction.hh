/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/validation/include/ValidationDetectorConstruction.hh
/// \brief Definition of the ValidationDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//
// 20260109  M. Kelsey -- G4CMP-569:  Remove unused electrodeSensitivity.

#ifndef ValidationDetectorConstruction_h
#define ValidationDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LatticeLogical.hh"
#include "G4CMPSurfaceProperty.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;


class ValidationDetectorConstruction : public G4VUserDetectorConstruction {
public:
  ValidationDetectorConstruction();
  virtual ~ValidationDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty* surfProp);

private:
  G4Material* fLiquidHelium;
  G4Material* fSilicon;
  G4Material* fAluminum;
  G4Material* fGermanium;
  G4Material* fTungsten;
  G4Material* fCopper;
  G4Material* fNiobium;
  G4VPhysicalVolume* fWorldPhys;

  //Interfaces for Geometry 1
  G4CMPSurfaceProperty* fVacSurfProp;
  G4CMPSurfaceProperty* fSiGeSurfProp;
  G4CMPSurfaceProperty* fGeAl1SurfProp;
  G4CMPSurfaceProperty* fAl1Al2FromAboveSurfProp;
  G4CMPSurfaceProperty* fAl1Al2FromBelowSurfProp;
  G4CMPSurfaceProperty* fAl2NbSurfProp;
  G4CMPSurfaceProperty* fAl2Al3SurfProp;
  G4CMPSurfaceProperty* fAl3NbSurfProp;

  G4bool fConstructed;		// Flag to not re-recreate surface properties

  //Interfaces for more complicated Example 7
  std::map<std::string,G4CMPSurfaceProperty*> fBorderContainer;
  std::map<std::string,G4LatticeLogical*> fLogicalLatticeContainer;
  G4CMPSurfaceProperty* fVacVacInterface;
  G4CMPSurfaceProperty* fAlVacInterface;
  G4CMPSurfaceProperty* fAlAlInterface;
  G4CMPSurfaceProperty* fSiAlInterface;
  G4CMPSurfaceProperty* fSiCuInterface;
  G4CMPSurfaceProperty* fSiVacInterface;
  G4CMPSurfaceProperty* fCuVacInterface;
};

#endif
