#ifndef Tst1PhysicsList_h
#define Tst1PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//#include "LatticeManager2.hh"

class Tst1PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst1PhysicsList();
   ~Tst1PhysicsList();

  public:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
private:
   
};

#endif



