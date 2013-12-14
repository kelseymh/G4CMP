
#ifndef Tst1PrimaryGeneratorAction_h
#define Tst1PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "LatticeManager2.hh"


class G4ParticleGun;
class G4Event;

class Tst1PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst1PrimaryGeneratorAction();    
  virtual ~Tst1PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                particleGun;

};


#endif


