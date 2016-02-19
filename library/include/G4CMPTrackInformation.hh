#ifndef G4CMPTrackInfo_h
#define G4CMPTrackInfo_h 1

#include "G4VAuxiliaryTrackInformation.hh"
#include "G4ThreeVector.hh"

class G4CMPTrackInformation : public G4VAuxiliaryTrackInformation {
public:
  // Default ctor probably shouldn't be used.
  G4CMPTrackInformation();

  // Phonon ctor
  G4CMPTrackInformation(G4ThreeVector k);

  // Charge ctor
  G4CMPTrackInformation(G4double l, G4double m, G4int v);

  void inline SetPhononK(G4ThreeVector k)                    { phononKVec = k; }
  void inline SetValleyIndex(G4int valley)         { chargeValleyIdx = valley; }
  void inline SetEffectiveMass(G4double m)                           { mc = m; }
  void inline SetScatterLength(G4double l)                           { l0 = l; }
  void inline SetReflectionCount(G4int n)                     { reflCount = n; }
  void inline IncrementReflectionCount()                        { ++reflCount; }
  G4ThreeVector inline GetPhononK() const                 { return phononKVec; }
  G4int inline GetValleyIndex() const                { return chargeValleyIdx; }
  G4double inline GetEffectiveMass() const                        { return mc; }
  G4double inline GetScatterLength() const                        { return l0; }
  G4int inline GetReflectionCount()                        { return reflCount; }

  virtual void Print() const override;

private:
  G4double l0;
  G4double mc;
  G4ThreeVector phononKVec;
  G4int chargeValleyIdx;
  G4int reflCount;
};

#endif
