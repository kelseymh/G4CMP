#include "G4CMPTrackInformation.hh"
#include "G4ThreeVector.hh"

// Default values are non-physical
G4CMPTrackInformation::G4CMPTrackInformation() :
  l0(0.), mc(0.), phononKVec(G4ThreeVector()), chargeValleyIdx(-1)
{}

G4CMPTrackInformation::G4CMPTrackInformation(G4ThreeVector k) :
  l0(0.), mc(0.), phononKVec(k), chargeValleyIdx(-1)
{}

G4CMPTrackInformation::G4CMPTrackInformation(G4double l,
                                             G4double m, G4int v) :
  l0(l), mc(m), phononKVec(G4ThreeVector()), chargeValleyIdx(v)
{}

void G4CMPTrackInformation::Print() const
{
//TODO
}
