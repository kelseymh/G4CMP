/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20160831  M. Kelsey -- Add optional electrode geometry class
// 20160907  M. Kelsey -- Protect against (allowed!) null electrode pointers
// 20170627  M. Kelsey -- Take ownership of electrode pointers and delete
// 20190806  M. Kelsey -- Add local data for frequency-dependent scattering
//		probabilities, and computation functions.

#include "G4CMPSurfaceProperty.hh"
#include "G4CMPVElectrodePattern.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <functional>
#include <vector>


// Constructors and destructor

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                                           G4SurfaceType stype)
  : G4SurfaceProperty(name, stype), theChargeElectrode(0),
    thePhononElectrode(0), anharmonicMaxFreq(0.), diffuseMaxFreq(0.) {;}

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

  // TEMPORARY: Hardcode phonon surface scattering parameters
  // Are these material dependent?  Are they device dependent?
  const G4double GHz = 1e9 * hertz;
  AddSurfaceAnharmonicCutoff(520*GHz);
  AddSurfaceDiffuseCutoff(350*GHz);
  AddSurfaceAnharmonicCoeffs({{0,0,0,0,0,1.51e-14}}, GHz);
  AddDiffuseReflectionCoeffs({{5.88e-2,7.83e-4,-2.47e-6,1.71e-8,-2.98e-11}}, GHz);
  AddSpecularReflectionCoeffs({{0.928,-2.03e-4,-3.21e-6,3.1e-9,2.9e-13}}, GHz);
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


// Frequency dependent phonon surface scattering probabilities

void G4CMPSurfaceProperty::
SaveCoeffs(std::vector<G4double>& buffer,
	   const std::vector<G4double>& coeff, G4double units) {
  buffer = coeff;
  if (units > 0.) {
    G4double unitpow = 1.;
    for (size_t i=0; i<buffer.size(); i++) {
      buffer[i] *= unitpow;		// Each coefficient gets units^i
      unitpow *= units;
    }
  }
}


// Compute phonon surface scattering probabilities

G4double G4CMPSurfaceProperty::
ExpandCoeffsPoly(G4double freq, const std::vector<G4double>& coeff) const {
  // Polynomial expansion like ((a[3]*x + a[2])*x + a[1])*x + a[0]
  G4double prob = 0.;
  for (size_t i=coeff.size(); i>0; prob=prob*freq+coeff[--i]);

  return prob;
}

G4double G4CMPSurfaceProperty::AnharmonicReflProb(G4double freq) const {
  if (freq > anharmonicMaxFreq) return 0.;

  return ExpandCoeffsPoly(freq, anharmonicCoeffs);
}

G4double G4CMPSurfaceProperty::DiffuseReflProb(G4double freq) const {
  if (freq > anharmonicMaxFreq)
    return 1. - thePhononMatPropTable.GetConstProperty("specProb");

  if (freq > diffuseMaxFreq) freq = diffuseMaxFreq;	// Flat response

  return ExpandCoeffsPoly(freq, diffuseCoeffs);
}

G4double G4CMPSurfaceProperty::SpecularReflProb(G4double freq) const {
  if (freq > diffuseMaxFreq)
    return 1. - DiffuseReflProb(freq) - AnharmonicReflProb(freq);

  return ExpandCoeffsPoly(freq, specularCoeffs);
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
