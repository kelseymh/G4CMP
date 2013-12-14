#ifndef LPhonon_h
#define LPhonon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"


class LPhonon : public G4ParticleDefinition
{
 private:
   static LPhonon* theInstance;

 private:
  LPhonon () {}


 public:

  
  /*  static void setLatticeManager(LatticeManager2*);
      static LatticeManager2* getLatticeManager();*/
  ~LPhonon (){}

  static LPhonon* Definition();
  static LPhonon* PhononDefinition();

};

#endif

