/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20160831  M. Kelsey -- Add optional electrode geometry class
// 20170525  M. Kelsey -- Add default "rule of five" copy/move operators
// 20170627  M. Kelsey -- Return non-const pointers, for functional use
// 20200601  G4CMP-206: Need thread-local copies of electrode pointers

#ifndef G4CMPSurfaceProperty_h
#define G4CMPSurfaceProperty_h 1

#include "G4SurfaceProperty.hh"
#include "G4MaterialPropertiesTable.hh"
#include <map>


class G4CMPVElectrodePattern;


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

  G4CMPSurfaceProperty(const G4CMPSurfaceProperty&) = default;
  G4CMPSurfaceProperty(G4CMPSurfaceProperty&&) = default;
  G4CMPSurfaceProperty& operator=(const G4CMPSurfaceProperty&) = default;
  G4CMPSurfaceProperty& operator=(G4CMPSurfaceProperty&&) = default;

  G4bool operator==(const G4SurfaceProperty &right) const;
  G4bool operator!=(const G4SurfaceProperty &right) const;

  // Accessors for charge-pair and phonon boundary parameters
  // NOTE:  Must return non-const pointer as Tables can't be const
  G4MaterialPropertiesTable* GetChargeMaterialPropertiesTablePointer() const {
    return const_cast<G4MaterialPropertiesTable*>(&theChargeMatPropTable);
  }

  G4MaterialPropertiesTable* GetPhononMaterialPropertiesTablePointer() const {
    return const_cast<G4MaterialPropertiesTable*>(&thePhononMatPropTable);
  }

  // NOTE:  These return by value because Tables can't be const
  G4MaterialPropertiesTable
  GetChargeMaterialPropertiesTable() const { return theChargeMatPropTable; }

  G4MaterialPropertiesTable
  GetPhononMaterialPropertiesTable() const { return thePhononMatPropTable; }

  // Accessors to fill charge-pair and phonon boundary parameters
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable& mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable& mpt);

  void FillChargeMaterialPropertiesTable(G4double qAbsProb, G4double qReflProb,
                                         G4double eMinK,    G4double hMinK);

  void FillPhononMaterialPropertiesTable(G4double pAbsProb,  G4double pReflProb,
                                         G4double pSpecProb, G4double pMinK);

  // Complex electrode geometries
  void SetChargeElectrode(G4CMPVElectrodePattern* cel);
  void SetPhononElectrode(G4CMPVElectrodePattern* pel);

  // Accessors, used by worker threads
  G4CMPVElectrodePattern* GetChargeElectrode() const;
  G4CMPVElectrodePattern* GetPhononElectrode() const;

  virtual void DumpInfo() const;	// To be implemented

protected:
  G4MaterialPropertiesTable theChargeMatPropTable;
  G4MaterialPropertiesTable thePhononMatPropTable;

  G4CMPVElectrodePattern* theChargeElectrode;
  G4CMPVElectrodePattern* thePhononElectrode;

  // These lists will be pre-allocated, with values entered by thread
  mutable std::map<G4int, G4CMPVElectrodePattern*> workerChargeElectrode;
  mutable std::map<G4int, G4CMPVElectrodePattern*> workerPhononElectrode;

  // These args should be const, but G4MaterialPropertiesTables is silly.
  G4bool IsValidChargePropTable(G4MaterialPropertiesTable& propTab) const;
  G4bool IsValidPhononPropTable(G4MaterialPropertiesTable& propTab) const;
};

#endif
