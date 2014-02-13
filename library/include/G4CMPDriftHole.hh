#ifndef G4CMPDriftHole_h
#define G4CMPDriftHole_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class G4CMPDriftHole : public G4ParticleDefinition{

private:
  static G4CMPDriftHole* theInstance;

private:
  G4CMPDriftHole(){}

public:
  ~G4CMPDriftHole(){}

  static G4CMPDriftHole* Definition();
  static G4CMPDriftHole* G4CMPDriftHoleDefinition();
};






#endif
