#ifndef Tst1PhysicsList_h
#define Tst1PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "ChannelingPhysicsMessenger.hh"

class ChannelingPhysicsList: public G4VUserPhysicsList
{
private:
    G4String fFileName;
    G4String fScatteringType;

    G4bool bWrapperOn;
    G4bool bChannelingOn;
    G4bool bDecayOn;

    ChannelingPhysicsMessenger *fMessenger;
    
public:
    ChannelingPhysicsList();
    ~ChannelingPhysicsList();
    
    //Add processes
    void AddStandardSS(G4ParticleDefinition*);
    void AddStandardNR(G4ParticleDefinition*);
    void AddChanneling();
    void AddDecay();
    void AddEM();
    void AddInelaticProcesses();
    void AddChannelingStepLimiter();
    
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
    
    void SetFileName(const G4String&);
    G4String GetFileName();

    void SetScatteringType(const G4String&);
    G4String GetScatteringType();

    void EnableChanneling(G4bool flag);
    G4bool GetChannelingState() {return bChannelingOn;};
    
    void EnableWrapper(G4bool flag);
    G4bool GetWrapperState() {return bWrapperOn;};
    
    void EnableDecay(G4bool flag);
    G4bool GetDecayState() {return bDecayOn;};
    
};

#endif
