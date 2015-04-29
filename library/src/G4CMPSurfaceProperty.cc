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
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4CMPSurfaceProperty.hh"

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
           G4SurfaceType type) : G4SurfaceProperty(name,type)
{;}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4String& name,
                         G4double prob, G4double deltaV,
                         G4double minKe, G4double minKh,
                         G4double V,
                         G4SurfaceType type)
                         : G4SurfaceProperty(name,type), absProb(prob),
                         absDeltaV(deltaV), minKElec(minKe), minKHole(minKh),
                         electrodeV(V)
{;}


G4CMPSurfaceProperty::~G4CMPSurfaceProperty()
{;}

G4CMPSurfaceProperty::G4CMPSurfaceProperty(const G4CMPSurfaceProperty &right)
  : G4SurfaceProperty(right.theName,right.theType)
{
    *this = right;
    this->theName    = right.theName;
    this->theType    = right.theType;
    this->absProb    = right.absProb;
    this->absDeltaV  = right.absDeltaV;
    this->minKElec   = right.minKElec;
    this->minKHole   = right.minKHole;
    this->electrodeV = right.electrodeV;
}

G4CMPSurfaceProperty& G4CMPSurfaceProperty::operator=(const G4CMPSurfaceProperty& right)
{
  if (this != &right) {
      theName    = right.theName;
      type       = right.type;
      absProb    = right.absProb;
      absDeltaV  = right.absDeltaV;
      minKElec   = right.minKElec;
      minKHole   = right.minKHole;
      electrodeV = right.electrodeV;
  }
  return *this;
}

G4int G4CMPSurfaceProperty::operator==(const G4CMPSurfaceProperty &right) const
{
  return (this == (G4CMPSurfaceProperty *) &right);
}

G4int G4CMPSurfaceProperty::operator!=(const G4CMPSurfaceProperty &right) const
{
  return (this != (G4CMPSurfaceProperty *) &right);
}
        ////////////
        // Methods
        ////////////

void G4CMPSurfaceProperty::DumpInfo() const
{

  // Dump info for surface
  // TO DO

}
