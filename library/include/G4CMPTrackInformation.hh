#ifndef G4CMPTrackInfo_h
#define G4CMPTrackInfo_h 1

#include "G4VAuxiliaryTrackInformation.hh"
#include "G4ThreeVector.hh"

class G4CMPTrackInformation : public G4VAuxiliaryTrackInformation {
public:
  G4CMPTrackInformation();

  void inline SetK(G4ThreeVector k)           {phononKVec = k;}
  void inline SetValleyIndex(G4int valley)    {chargeValleyIdx = valley;}
  //void inline SetReflectionCount(G4int n)     {chargeReflCount = n;}
  //void inline IncrementReflectionCount()      {++chargeReflCount;}
  G4ThreeVector inline GetK()                 {return phononKVec;}
  G4int inline GetValleyIndex()               {return chargeValleyIdx;}

  virtual void Print() const override;

private:
  G4ThreeVector phononKVec;
  G4int chargeValleyIdx;
  //G4int chargeReflCount;
};

#endif
