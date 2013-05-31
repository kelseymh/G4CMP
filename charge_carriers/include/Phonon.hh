
#ifndef Phonon_h
#define Phonon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class Phonon : public G4ParticleDefinition
{
 private:
   static Phonon* theInstance;

 
  Phonon () {}

  public:
   ~Phonon (){}

   static Phonon* Definition();
   static Phonon* PhononDefinition();

};

#endif

