/* Header File for AnharmonicDecay utility class */

#ifndef G4CMPAnharmonicDecay_h
#define G4CMPAnharmonicDecay_h

#include "G4CMPProcessUtils.hh"

class G4CMPAnharmonicDecay : public G4CMPProcessUtils {
public:
  G4CMPAnharmonicDecay();
  virtual ~G4CMPAnharmonicDecay();
  virtual void DoDecay(const G4Track&, const G4Step&,
                                    G4ParticleChange*);

private:
  inline double GetLTDecayProb(G4double, G4double) const;
  inline double GetTTDecayProb(G4double, G4double) const;
  inline double MakeLDeviation(G4double, G4double) const;
  inline double MakeTTDeviation(G4double, G4double) const;
  inline double MakeTDeviation(G4double, G4double) const;

  void MakeTTSecondaries(const G4Track&, G4ParticleChange);
  void MakeLTSecondaries(const G4Track&, G4ParticleChange);

  G4double fBeta, fGamma, fLambda, fMu; // Local buffers for calculations
  G4double fvLvT; // Ratio of sound speeds

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif /* G4CMPAnharmonicDecay_h */
