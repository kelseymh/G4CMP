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

#include "G4ios.hh"
#include "globals.hh"
#include "G4CMPSurfaceProperty.hh"

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
           G4SurfaceType stype) : G4SurfaceProperty(name, stype)
{;}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                         G4double qAbsprob, G4double V,
                         G4double deltaV, G4double minKe,
                         G4double minKh, G4double pAbsprob, G4double specProb,
                         G4double gapEnergy, G4double lowQPLimit,
                         G4double pLifetime, G4double vSound,
                         G4double lifetimeVsESlope, G4double filmThickness,
                         G4SurfaceType stype)
                         : G4SurfaceProperty(name, stype)
{
  theChargeMatPropTable.AddConstProperty("absProb", qAbsprob);
  theChargeMatPropTable.AddConstProperty("electrodeV", V);
  theChargeMatPropTable.AddConstProperty("absDeltaV", deltaV);
  theChargeMatPropTable.AddConstProperty("minKElec", minKe);
  theChargeMatPropTable.AddConstProperty("minKHole", minKh);

  thePhononMatPropTable.AddConstProperty("absProb", pAbsprob);
  thePhononMatPropTable.AddConstProperty("specProb", specProb);
  thePhononMatPropTable.AddConstProperty("gapEnergy", gapEnergy);
  thePhononMatPropTable.AddConstProperty("lowQPLimit", lowQPLimit);
  thePhononMatPropTable.AddConstProperty("phononLifetime", pLifetime);
  thePhononMatPropTable.AddConstProperty("phononLifetimeSlope", lifetimeVsESlope);
  thePhononMatPropTable.AddConstProperty("vSound", vSound);
  thePhononMatPropTable.AddConstProperty("filmThickness", filmThickness);
}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                         G4double qAbsprob, G4double V,
                         G4double deltaV, G4double minKe,
                         G4double minKh,
                         G4SurfaceType stype)
                         : G4SurfaceProperty(name,stype)
{
  theChargeMatPropTable.AddConstProperty("absProb", qAbsprob);
  theChargeMatPropTable.AddConstProperty("electrodeV", V);
  theChargeMatPropTable.AddConstProperty("absDeltaV", deltaV);
  theChargeMatPropTable.AddConstProperty("minKElec", minKe);
  theChargeMatPropTable.AddConstProperty("minKHole", minKh);
}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                         G4double pAbsprob, G4double specProb,
                         G4double gapEnergy, G4double lowQPLimit,
                         G4double pLifetime, G4double lifetimeVsESlope,
                         G4double vSound, G4double filmThickness,
                         G4SurfaceType stype)
                         : G4SurfaceProperty(name, stype)
{
  thePhononMatPropTable.AddConstProperty("absProb", pAbsprob);
  thePhononMatPropTable.AddConstProperty("specProb", specProb);
  thePhononMatPropTable.AddConstProperty("gapEnergy", gapEnergy);
  thePhononMatPropTable.AddConstProperty("lowQPLimit", lowQPLimit);
  thePhononMatPropTable.AddConstProperty("phononLifetime", pLifetime);
  thePhononMatPropTable.AddConstProperty("phononLifetimeSlope", lifetimeVsESlope);
  thePhononMatPropTable.AddConstProperty("vSound", vSound);
  thePhononMatPropTable.AddConstProperty("filmThickness", filmThickness);
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
