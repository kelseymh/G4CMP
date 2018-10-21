/* Header File for AnharmonicDecay utility class */

#ifndef G4CMPAnharmonicDecay_h
#define G4CMPAnharmonicDecay_h

#include "G4CMPProcessUtils.hh"
#include <iosfwd>

class G4ParticleChange;
class G4Step;
class G4Track;


class G4CMPAnharmonicDecay : public G4CMPProcessUtils {
public:
  G4CMPAnharmonicDecay(G4int vb=0);
  virtual ~G4CMPAnharmonicDecay() {;}

  void SetVerboseLevel(G4int vb) { verboseLevel = vb;}
  G4int GetVerboseLevel() const { return verboseLevel; }

  void DoDecay(const G4Track&, const G4Step&, G4ParticleChange&);

private:
  G4double GetLTDecayProb(G4double, G4double) const;
  G4double GetTTDecayProb(G4double, G4double) const;
  G4double MakeLDeviation(G4double, G4double) const;
  G4double MakeTTDeviation(G4double, G4double) const;
  G4double MakeTDeviation(G4double, G4double) const;

  void MakeTTSecondaries(const G4Track&, G4ParticleChange&);
  void MakeLTSecondaries(const G4Track&, G4ParticleChange&);

  G4int verboseLevel;			// For diagnostic output
  G4double fBeta, fGamma, fLambda, fMu; // Local buffers for decay parameters
  G4double fvLvT; 			// Ratio of sound speeds

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif /* G4CMPAnharmonicDecay_h */
