#ifndef G4CMPDriftElectron_h
#define G4CMPDriftElectron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class G4CMPDriftElectron : public G4ParticleDefinition{

private:
  static G4CMPDriftElectron* theInstance;

private:
  G4CMPDriftElectron(){}

public:
  ~G4CMPDriftElectron(){}

  static G4CMPDriftElectron* Definition();
  static G4CMPDriftElectron* G4CMPDriftElectronDefinition();
};






#endif
