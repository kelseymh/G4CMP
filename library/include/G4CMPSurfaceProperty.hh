//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$

#ifndef G4CMPSurfaceProperty_h
#define G4CMPSurfaceProperty_h 1

#include "G4Types.hh"
#include "G4Physics2DVector.hh"
#include "G4SurfaceProperty.hh"
#include "G4MaterialPropertiesTable.hh"

class G4CMPSurfaceProperty : public G4SurfaceProperty {
public:
  //Empty constructor. Not very useful.
  G4CMPSurfaceProperty(const G4String& name,
                         G4SurfaceType stype = dielectric_dielectric);

  //Full constructor
  G4CMPSurfaceProperty(const G4String& name,
                       G4double qAbsProb, //Prob. to absorb charge no matter what
                       G4double V, //Voltage of electrode
                       G4double deltaV, //Absorb charge if voltage >= |V-deltaV|
                       G4double minKe, //Min wave number to absorb electron
                       G4double minKh, //Min wave number to absorb hole
                       G4double pAbsProb, //Prob. to absorb phonon
                       G4double specProb, //Prob. of specular reflection
                       G4double gapEnergy, //Band gap energy of second surf.
                       G4double lowQPLimit, //Down convert phonon until energy < lowQPLimit*gapEnergy
                       G4double pLifetime, //Average phonon lifetime in file @ E=2*gapE
                       G4double vSound, //Speed of sound in film
                       G4double lifetimeVsESlope, //Unitless slope to define tau(E)
                       G4double filmThickness,
                       G4SurfaceType stype = dielectric_dielectric);

  //Charge-only constructor for convenience
  G4CMPSurfaceProperty(const G4String& name,
                       G4double qAbsprob, G4double V,
                       G4double deltaV, G4double minKe,
                       G4double minKh,
                       G4SurfaceType stype = dielectric_dielectric);

  //Phonon-only consutrctor for convenience
  G4CMPSurfaceProperty(const G4String& name,
                       G4double pAbsprob, G4double specProb,
                       G4double gapEnergy,G4double lowQPLimit,
                       G4double pLifetime, G4double vSound,
                       G4double lifetimeVsESlope, G4double filmThickness,
                       G4SurfaceType stype = dielectric_dielectric);

  G4int operator==(const G4CMPSurfaceProperty &right) const;
  G4int operator!=(const G4CMPSurfaceProperty &right) const;

  inline const G4MaterialPropertiesTable* GetChargeMaterialPropertiesTablePointer() const
                       { return &theChargeMatPropTable; }
  inline const G4MaterialPropertiesTable* GetPhononMaterialPropertiesTablePointer() const
                       { return &thePhononMatPropTable; }
  inline G4MaterialPropertiesTable GetChargeMaterialPropertiesTable() const
                       { return theChargeMatPropTable; }
  inline G4MaterialPropertiesTable GetPhononMaterialPropertiesTable() const
                       { return thePhononMatPropTable; }

  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable mpt);

  void DumpInfo() const;	// To be implemented

private:
  G4MaterialPropertiesTable theChargeMatPropTable;
  G4MaterialPropertiesTable thePhononMatPropTable;

  G4MaterialPropertiesTable CopyMaterialPropertiesTable(
    G4MaterialPropertiesTable *mpt);
};

#endif
