/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/quasiparticle/include/QuasiparticleDetectorConstruction.hh
/// \brief Definition of the QuasiparticleDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//
// 20260109  M. Kelsey -- G4CMP-569:  Remove unused electrodeSensitivity.

#ifndef QuasiparticleDetectorConstruction_h
#define QuasiparticleDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LatticeLogical.hh"
#include "G4CMPSurfaceProperty.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;


class QuasiparticleDetectorConstruction : public G4VUserDetectorConstruction {
public:
  QuasiparticleDetectorConstruction();
  virtual ~QuasiparticleDetectorConstruction();
  
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
  G4Material* fTungsten;
  G4Material* fCopper;
  G4Material* fNiobium;
  G4VPhysicalVolume* fWorldPhys;

  //Interfaces for Examples 1-6
  G4CMPSurfaceProperty* topSurfProp; //Al/Si
  G4CMPSurfaceProperty* topSurfProp2; //Nb/Si
  G4CMPSurfaceProperty* vacSurfProp;
  G4CMPSurfaceProperty* wallSurfProp;
  G4CMPSurfaceProperty* copTopSurfProp; //Al/Cu
  G4CMPSurfaceProperty* copTopSurfProp2; //Nb/Cu
  G4CMPSurfaceProperty* alNbSurfProp; //Al/Nb

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

