#ifndef DriftingHole_h
#define DriftingHole_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class DriftingHole : public G4ParticleDefinition{

private:
  static DriftingHole* theInstance;

private:
  DriftingHole(){}

public:
  ~DriftingHole(){}

  static DriftingHole* Definition();
  static DriftingHole* DriftingHoleDefinition();
};






#endif
