/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20160831  M. Kelsey -- Add optional electrode geometry class
// 20160907  M. Kelsey -- Protect against (allowed!) null electrode pointers
// 20170627  M. Kelsey -- Take ownership of electrode pointers and delete

#include "G4CMPSurfaceProperty.hh"
#include "G4CMPVElectrodePattern.hh"


// Constructors and destructor

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                                           G4SurfaceType stype)
  : G4SurfaceProperty(name, stype), theChargeElectrode(0),
    thePhononElectrode(0) {;}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                                           G4double qAbsProb,
                                           G4double qReflProb,
                                           G4double eMinK,
                                           G4double hMinK,
                                           G4double pAbsProb,
                                           G4double pReflProb,
                                           G4double pSpecProb,
                                           G4double pMinK,
                                           G4SurfaceType stype)
: G4CMPSurfaceProperty(name, stype) {
  FillChargeMaterialPropertiesTable(qAbsProb, qReflProb, eMinK, hMinK);
  FillPhononMaterialPropertiesTable(pAbsProb, pReflProb, pSpecProb, pMinK);
}

G4CMPSurfaceProperty::~G4CMPSurfaceProperty() {
  // Find ourselves in the static SurfaceProperty table and remove ourselves.
  std::vector<G4SurfaceProperty*>::iterator me =
    std::find_if(theSurfacePropertyTable.begin(), theSurfacePropertyTable.end(),
                 [this](G4SurfaceProperty* a)->G4bool {return *this == *a;});
  if (me != theSurfacePropertyTable.end())
    theSurfacePropertyTable.erase(me);

  // Delete electrodes associated with this surface
  delete theChargeElectrode; theChargeElectrode=0;
  delete thePhononElectrode; thePhononElectrode=0;
}

G4bool G4CMPSurfaceProperty::operator==(const G4SurfaceProperty& right) const {
  return (this == dynamic_cast<const G4CMPSurfaceProperty*>(&right));
}

G4bool G4CMPSurfaceProperty::operator!=(const G4SurfaceProperty &right) const {
  return !(*this == right);
}

void G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable(
                            G4MaterialPropertiesTable* mpt) {
  if (IsValidChargePropTable(*mpt)) {
    theChargeMatPropTable = *mpt;
  } else {
    G4Exception("G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable",
                "detector001", RunMustBeAborted,
                "Tried to set properties table to one that is not valid.");
  }
}

void G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable(
                            G4MaterialPropertiesTable* mpt) {
  if (IsValidChargePropTable(*mpt)) {
    thePhononMatPropTable = *mpt;
  } else {
    G4Exception("G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable",
                "detector002", RunMustBeAborted,
                "Tried to set properties table to one that is not valid.");
  }
}

void G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable(
  G4MaterialPropertiesTable& mpt) {
  if (IsValidChargePropTable(mpt)) {
    theChargeMatPropTable = mpt;
  } else {
    G4Exception("G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable",
                "detector003", RunMustBeAborted,
                "Tried to set properties table to one that is not valid.");
  }
}

void G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable(
  G4MaterialPropertiesTable& mpt) {
  if (IsValidChargePropTable(mpt)) {
    thePhononMatPropTable = mpt;
  } else {
    G4Exception("G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable",
                "detector004", RunMustBeAborted,
                "Tried to set properties table to one that is not valid.");
  }
}

void G4CMPSurfaceProperty::FillChargeMaterialPropertiesTable(G4double qAbsProb,
                                                             G4double qReflProb,
                                                             G4double eMinK,
                                                             G4double hMinK) {
  theChargeMatPropTable.AddConstProperty("absProb", qAbsProb);
  theChargeMatPropTable.AddConstProperty("reflProb", qReflProb);
  theChargeMatPropTable.AddConstProperty("minKElec", eMinK);
  theChargeMatPropTable.AddConstProperty("minKHole", hMinK);
}

void G4CMPSurfaceProperty::FillPhononMaterialPropertiesTable(G4double pAbsProb,
                                                             G4double pReflProb,
                                                             G4double pSpecProb,
                                                             G4double pMinK) {
  thePhononMatPropTable.AddConstProperty("absProb", pAbsProb);
  thePhononMatPropTable.AddConstProperty("reflProb", pReflProb);
  thePhononMatPropTable.AddConstProperty("specProb", pSpecProb);
  thePhononMatPropTable.AddConstProperty("absMinK", pMinK);
}


// Complex electrode geometries

void G4CMPSurfaceProperty::SetChargeElectrode(G4CMPVElectrodePattern* cel) {
  theChargeElectrode = cel;
  if (cel) theChargeElectrode->UseSurfaceTable(&theChargeMatPropTable);
}

void G4CMPSurfaceProperty::SetPhononElectrode(G4CMPVElectrodePattern* pel) {
  thePhononElectrode = pel;
  if (pel) thePhononElectrode->UseSurfaceTable(&thePhononMatPropTable);
}


// Report configuration parameters for diagnostics

void G4CMPSurfaceProperty::DumpInfo() const {
  // FIXME:  Stupid Tables don't have any const accessors!
  const_cast<G4MaterialPropertiesTable*>(&theChargeMatPropTable)->DumpTable();
  const_cast<G4MaterialPropertiesTable*>(&thePhononMatPropTable)->DumpTable();
}

G4bool G4CMPSurfaceProperty::
IsValidChargePropTable(G4MaterialPropertiesTable& propTab) const {
  // A property table is valid for us if it at least contains all of the
  // properties that we require.
  return (propTab.ConstPropertyExists("absProb") &&
          propTab.ConstPropertyExists("reflProb") &&
          propTab.ConstPropertyExists("minKElec") &&
          propTab.ConstPropertyExists("minKHole"));
}

G4bool G4CMPSurfaceProperty::
IsValidPhononPropTable(G4MaterialPropertiesTable& propTab) const {
  // A property table is valid for us if it at least contains all of the
  // properties that we require.
  return (propTab.ConstPropertyExists("absProb") &&
          propTab.ConstPropertyExists("reflProb") &&
          propTab.ConstPropertyExists("specProb") &&
          propTab.ConstPropertyExists("absMinK"));
}
