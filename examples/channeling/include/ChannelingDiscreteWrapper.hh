/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChannelingDiscreteWrapper_h
#define ChannelingDiscreteWrapper_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "ChannelingParticleUserInfo.hh"

class G4Material;

class ChannelingDiscreteWrapper : public G4VDiscreteProcess
{
public:
    
    ChannelingDiscreteWrapper(const G4String& processName ="ChannelingDiscreteWrapper" );
    ChannelingDiscreteWrapper(const G4String& , G4VDiscreteProcess*);
    
    virtual ~ChannelingDiscreteWrapper();
    
public:
    void RegisterProcess(G4VDiscreteProcess*);
    void RegisterProcess(G4VDiscreteProcess*,G4int);
    
    G4double GetDensity(const G4Track&);
    G4double GetDensityPreviousStep(const G4Track&);
    
    void SetNucleiOrElectronFlag(G4int);
    G4int GetNucleiOrElectronFlag();
    
private:
    // hide assignment operator as private
    ChannelingDiscreteWrapper(ChannelingDiscreteWrapper&);
    ChannelingDiscreteWrapper& operator=(const ChannelingDiscreteWrapper&);
    
    //private data members
    G4int bNucleiOrElectronFlag; //Decide whether to use nuclei (+1) or electron (-1) or both (0) density to change parameters
    G4VDiscreteProcess* fRegisteredProcess;

    
    /////////////////////////////////////////////////////////
    /////////////////// GEANT4 PROCESS METHODS //////////////
    /////////////////////////////////////////////////////////
public:
    // DO IT
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
    
    // GPIL
    virtual G4double PostStepGetPhysicalInteractionLength (const G4Track&, G4double, G4ForceCondition*);
    
    // GENERAL
    void StartTracking(G4Track*);
    virtual G4bool IsApplicable(const G4ParticleDefinition&);

    // PHYSICS TABLE
    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    virtual void PreparePhysicsTable(const G4ParticleDefinition&);
    virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,const G4String&, G4bool);
    virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,const G4String&, G4bool);

protected:
    // MFP
    virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* );
 
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

};

#endif










