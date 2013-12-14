
#ifndef TPhononFast_h
#define TPhononFast_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "LatticeManager2.hh"


class TPhononFast : public G4ParticleDefinition
{
 private:
   static TPhononFast* theInstance;

 private:
  TPhononFast () {}
 

 public:
   ~TPhononFast (){}

  static LatticeManager2* LM;

  static void setLatticeManager(LatticeManager2*);
  static LatticeManager2* getLatticeManager();

  static TPhononFast* Definition();
  static TPhononFast* PhononDefinition();

};

#endif

