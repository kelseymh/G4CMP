/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$

#include "G4ios.hh"
#include "globals.hh"
#include "G4CMPSurfaceProperty.hh"

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                                           G4SurfaceType stype)
                                                : G4SurfaceProperty(name, stype)
{;}

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
                                                : G4SurfaceProperty(name, stype)
{
  theChargeMatPropTable.AddConstProperty("absProb", qAbsProb);
  theChargeMatPropTable.AddConstProperty("reflProb", qReflProb);
  theChargeMatPropTable.AddConstProperty("minKElec", eMinK);
  theChargeMatPropTable.AddConstProperty("minKHole", hMinK);

  thePhononMatPropTable.AddConstProperty("absProb", pAbsProb);
  thePhononMatPropTable.AddConstProperty("reflProb", pReflProb);
  thePhononMatPropTable.AddConstProperty("specProb", pSpecProb);
  thePhononMatPropTable.AddConstProperty("minK", pMinK);
}

G4int G4CMPSurfaceProperty::operator==(const G4CMPSurfaceProperty &right) const
{
  return (this == (G4CMPSurfaceProperty *) &right);
}

G4int G4CMPSurfaceProperty::operator!=(const G4CMPSurfaceProperty &right) const
{
  return (this != (G4CMPSurfaceProperty *) &right);
}

void G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable(
                            G4MaterialPropertiesTable* mpt)
{
  theChargeMatPropTable = CopyMaterialPropertiesTable(mpt);
}

void G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable(
                            G4MaterialPropertiesTable* mpt)
{
  thePhononMatPropTable = CopyMaterialPropertiesTable(mpt);
}

void G4CMPSurfaceProperty::SetChargeMaterialPropertiesTable(
  G4MaterialPropertiesTable mpt)
{
  theChargeMatPropTable = mpt;
}

void G4CMPSurfaceProperty::SetPhononMaterialPropertiesTable(
  G4MaterialPropertiesTable mpt)
{
  thePhononMatPropTable = mpt;
}

void G4CMPSurfaceProperty::FillChargeMaterialPropertiesTable(G4double qAbsProb,
                                                             G4double qReflProb,
                                                             G4double eMinK,
                                                             G4double hMinK)
{
  theChargeMatPropTable.AddConstProperty("absProb", qAbsProb);
  theChargeMatPropTable.AddConstProperty("reflProb", qReflProb);
  theChargeMatPropTable.AddConstProperty("minKElec", eMinK);
  theChargeMatPropTable.AddConstProperty("minKHole", hMinK);
}

void G4CMPSurfaceProperty::FillPhononMaterialPropertiesTable(G4double pAbsProb,
                                                             G4double pReflProb,
                                                             G4double pSpecProb,
                                                             G4double pMinK)
{
  thePhononMatPropTable.AddConstProperty("absProb", pAbsProb);
  thePhononMatPropTable.AddConstProperty("reflProb", pReflProb);
  thePhononMatPropTable.AddConstProperty("specProb", pSpecProb);
  thePhononMatPropTable.AddConstProperty("minK", pMinK);
}

void G4CMPSurfaceProperty::DumpInfo() const
{

  // Dump info for surface
  // TO DO

}

G4MaterialPropertiesTable G4CMPSurfaceProperty::CopyMaterialPropertiesTable(
  G4MaterialPropertiesTable* in)
{
  G4MaterialPropertiesTable out;
  const std::map<G4String,G4MaterialPropertyVector*,std::less<G4String> >*
    inMap = in->GetPropertiesMap();
  const std::map<G4String,G4double,std::less<G4String> >* inCMap =
                                              in->GetPropertiesCMap();
  std::map<G4String,G4MaterialPropertyVector*,std::less<G4String> >::const_iterator
   inItr = inMap->begin();
  std::map<G4String,G4double,std::less<G4String> >::const_iterator inCItr =
                                              inCMap->begin();

  for(; inItr != inMap->end(); ++inItr) {
    out.AddProperty(inItr->first, inItr->second);
  }

  for(; inCItr != inCMap->end(); ++inCItr) {
    out.AddConstProperty(inItr->first, inCItr->second);
  }

  return out;
}
