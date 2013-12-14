#ifndef DriftingElectron_h
#define DriftingElectron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class DriftingElectron : public G4ParticleDefinition{

private:
  static DriftingElectron* theInstance;

private:
  DriftingElectron(){}

public:
  ~DriftingElectron(){}

  static DriftingElectron* Definition();
  static DriftingElectron* DriftingElectronDefinition();
};






#endif
