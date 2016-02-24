/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 7649d43d17165e405358d4db730852c219d62318 $
//
//

#ifndef ChannelingUserInfo_h
#define ChannelingUserInfo_h 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

class ChannelingParticleUserInfo : public G4VUserTrackInformation
{
    
public:
    
    ChannelingParticleUserInfo();
    ~ChannelingParticleUserInfo();
    
    void SetCoherentEffect(G4int flag); 
    G4int HasBeenUnderCoherentEffect();
    
    void SetNucleiDensity(G4double);
    G4double GetNucleiDensity();
    
    void SetElectronDensity(G4double);
    G4double GetElectronDensity();
    
    G4double GetNucleiDensityPreviousStep();
    G4double GetElectronDensityPreviousStep();
    void StoreDensityPreviousStep();

    G4ThreeVector GetMomentumChanneled();
    void SetMomentumChanneled(G4ThreeVector);

    G4ThreeVector GetPositionChanneled();
    void SetPositionChanneled(G4ThreeVector);

    G4double GetEnergyChanneled();
    void SetEnergyChanneled(G4double);

    G4ThreeVector GetMomentumChanneledInitial();
    void SetMomentumChanneledInitial(G4ThreeVector);
    
    G4ThreeVector GetPositionChanneledInitial();
    void SetPositionChanneledInitial(G4ThreeVector);

    G4int GetNumberOfDechanneling();
    void IncreaseNumberOfDechanneling();
    

private:
    
    G4int bHasBeenUnderCoherentEffect; //Has been in channeling in the last step

    G4double fNucleiDensity; //Last value of density seen by channeled particle
    G4double fNucleiDensityPreviousStep;
    
    G4double fElectronDensity; //Last value of density seen by channeled particle
    G4double fElectronDensityPreviousStep;

    G4ThreeVector fMomentumInChanneling; //Last position of the particle in the channel
    G4ThreeVector fMomentumInChannelingInitial; //Last position of the particle in the channel

    G4ThreeVector fPositionInChanneling; //Last projection fof the particle momentum in the crystal reference system
    G4ThreeVector fPositionInChannelingInitial; //Last projection fof the particle momentum in the crystal reference system
    
    G4int fNumberOfDechanneling;
    

};



#endif
