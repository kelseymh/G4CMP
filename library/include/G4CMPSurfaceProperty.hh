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

#ifndef G4CMPSurfaceProperty_h
#define G4CMPSurfaceProperty_h 1

#include "G4Types.hh"
#include "G4Physics2DVector.hh"
#include "G4SurfaceProperty.hh"

// TODO: Someday will probably have to use MaterialPropertiesTable
//class G4MaterialPropertiesTable;

class G4CMPSurfaceProperty : public G4SurfaceProperty
{

public:
    G4CMPSurfaceProperty(const G4String& name,
                         G4SurfaceType type = dielectric_dielectric);

    G4CMPSurfaceProperty(const G4String& name,
                         G4double prob, G4double deltaV,
                         G4double minKe, G4double minKh,
                         G4double V, G4SurfaceType type = dielectric_dielectric);

    G4CMPSurfaceProperty(const G4CMPSurfaceProperty &right);
    G4CMPSurfaceProperty & operator=(const G4CMPSurfaceProperty &right);

    G4int operator==(const G4CMPSurfaceProperty &right) const;
    G4int operator!=(const G4CMPSurfaceProperty &right) const;

    virtual ~G4CMPSurfaceProperty();

    //inline G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
    //                   { return theMaterialPropertiesTable; }

    //inline void SetMaterialPropertiesTable(G4MaterialPropertiesTable *anMPT)
    //                   { theMaterialPropertiesTable = anMPT; }

    void DumpInfo() const;

    void SetType(const G4SurfaceType& type);
    inline G4SurfaceType GetType() {return type; }
    inline void SetAbsProb(G4double p) { absProb = p; }
    inline G4double GetAbsProb() { return absProb; }
    inline void SetAbsDeltaV(G4double p) { absDeltaV = p; }
    inline G4double GetAbsDeltaV() { return absDeltaV; }
    inline void SetMinKElec(G4double k) { minKElec = k; }
    inline G4double GetMinKElec() { return minKElec; }
    inline void SetMinKHole(G4double k) { minKHole = k; }
    inline G4double GetMinKHole() { return minKHole; }
    inline void SetElectrodeV(G4double v) { electrodeV = v; }
    inline G4double GetElectrodeV() { return electrodeV; }


    //void ReadLUTFile(void);

private:
    //G4MaterialPropertiesTable* theMaterialPropertiesTable;
    //FILE* readLUTFileHandle;
    G4SurfaceType type;
    G4double absProb;
    G4double absDeltaV;
    G4double minKElec;
    G4double minKHole;
    G4double electrodeV;
};

#endif /* G4OpticalSurface_h */
