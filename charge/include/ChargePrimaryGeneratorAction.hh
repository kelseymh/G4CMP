
#ifndef ChargePrimaryGeneratorAction_h
#define ChargePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4LatticeManager.hh"


class G4ParticleGun;
class G4Event;

class ChargePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ChargePrimaryGeneratorAction();    
  virtual ~ChargePrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                particleGun;

};


#endif


