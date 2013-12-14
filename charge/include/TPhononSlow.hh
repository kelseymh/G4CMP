
#ifndef TPhononSlow_h
#define TPhononSlow_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "LatticeManager2.hh"

class TPhononSlow : public G4ParticleDefinition
{
 private:
   static TPhononSlow* theInstance;

 private:
  TPhononSlow () {}

 public:

  static LatticeManager2* LM;
  
  static void setLatticeManager(LatticeManager2*);
  static LatticeManager2* getLatticeManager(); 
  ~TPhononSlow (){}

   static TPhononSlow* Definition();
   static TPhononSlow* PhononDefinition();

};

#endif

