/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$

#ifndef G4CMPSurfaceProperty_h
#define G4CMPSurfaceProperty_h 1

#include "G4SurfaceProperty.hh"
#include "G4MaterialPropertiesTable.hh"

class G4CMPSurfaceProperty : public G4SurfaceProperty {
public:
  // Empty constructor. Users must call at least one of the FillPropertiesTable
  // member functions. But, really, you shouldn't use this. It's dangerous and
  // I don't know why I put it here at all.
  G4CMPSurfaceProperty(const G4String& name,
                       G4SurfaceType stype = dielectric_dielectric);

  //Full constructor
  G4CMPSurfaceProperty(const G4String& name,
                       G4double qAbsProb, // Prob. to absorb charge carrier
                       G4double qReflProb, // If not absorbed, prob to reflect
                       G4double eMinK, //Min wave number to absorb electron
                       G4double hMinK, //Min wave number to absorb hole
                       G4double pAbsProb, // Prob. to absorb phonon
                       G4double pReflProb, // If not absorbed, prob to reflect
                       G4double pSpecProb, //Prob. of specular reflection
                       G4double pMinK, //Min wave number to absorb phonon
                       G4SurfaceType stype = dielectric_dielectric);

  virtual ~G4CMPSurfaceProperty();

  G4bool operator==(const G4SurfaceProperty &right) const;
  G4bool operator!=(const G4SurfaceProperty &right) const;

  const G4MaterialPropertiesTable* GetChargeMaterialPropertiesTablePointer() const
                       { return &theChargeMatPropTable; }
  const G4MaterialPropertiesTable* GetPhononMaterialPropertiesTablePointer() const
                       { return &thePhononMatPropTable; }
  G4MaterialPropertiesTable GetChargeMaterialPropertiesTable() const
                       { return theChargeMatPropTable; }
  G4MaterialPropertiesTable GetPhononMaterialPropertiesTable() const
                       { return thePhononMatPropTable; }

  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable mpt);

  void FillChargeMaterialPropertiesTable(G4double qAbsProb,
                                         G4double qReflProb,
                                         G4double eMinK,
                                         G4double hMinK);

  void FillPhononMaterialPropertiesTable(G4double pAbsProb,
                                         G4double pReflProb,
                                         G4double pSpecProb,
                                         G4double pMinK);

  void DumpInfo() const;	// To be implemented

protected:
  G4MaterialPropertiesTable theChargeMatPropTable;
  G4MaterialPropertiesTable thePhononMatPropTable;

  // These args should be const, but G4MaterialPropertiesTables is silly.
  virtual G4bool IsValidChargePropTable(G4MaterialPropertiesTable& propTab) const;
  virtual G4bool IsValidPhononPropTable(G4MaterialPropertiesTable& propTab) const;
};

#endif
