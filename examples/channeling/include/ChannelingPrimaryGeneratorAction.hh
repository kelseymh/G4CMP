/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: b119b8c1181ee0682e3666763a6a8dd0b21a8bc6 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ChannelingPrimaryGeneratorAction_h
#define ChannelingPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "ChannelingPrimaryGeneratorMessenger.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class ChannelingDetectorConstruction;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChannelingPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    ChannelingPrimaryGeneratorAction();
    virtual ~ChannelingPrimaryGeneratorAction();
    
    void GeneratePrimaries(G4Event*);

    void SetBeamDivergencDistribution(G4String div) {fDivDistribution = div;};
    G4String GetBeamDivergencDistribution() {return fDivDistribution;};

    void SetBeamDivergenceX(G4double div) {fDivX = div;};
    G4double GetBeamDivergenceX() {return fDivX;};

    void SetBeamDivergenceY(G4double div) {fDivY=div;};
    G4double GetBeamDivergenceY() {return fDivY;};

    void SetBeamCutX(G4double div) {fCutX = div;};
    G4double GetBeamCutX() {return fCutX;};
    
    void SetBeamCutY(G4double div) {fCutY = div;};
    G4double GetBeamCutY() {return fCutY;};

    void SetBeamWidthX(G4double width) {fWidthX = width;};
    G4double GetBeamWidthX() {return fWidthX;};
    
    void SetBeamWidthY(G4double width) {fWidthY=width;};
    G4double GetBeamWidthY() {return fWidthY;};

private:
    G4ParticleGun* fParticleGun;	 //pointer a to G4  class
    
    ChannelingPrimaryGeneratorMessenger *fMessenger;

    
    G4String fDivDistribution;
    G4double fCutX;
    G4double fCutY;
    G4double fDivX;
    G4double fDivY;
    G4double fWidthX;
    G4double fWidthY;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


