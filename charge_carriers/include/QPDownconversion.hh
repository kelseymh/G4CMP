#ifndef QPDownconversion_h
#define QPDownconversion_h 1

#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"

#include <iostream>
#include <fstream>

using namespace std;

class QPDownconversion{

public:

  QPDownconversion();
  ~QPDownconversion();

  void downconvert(const G4Track& aTrack, G4VParticleChange* aParticleChange, G4ThreeVector direction);

private:

  ofstream writer;

  G4double QPEnergy(G4double phononE);
  G4double PhEnergy(G4double QPEi);

};


#endif
