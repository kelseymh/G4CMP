#ifndef G4CMPInterValleyScattering_h
#define G4CMPInterValleyScattering_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class G4CMPInterValleyScattering : public G4VDiscreteProcess,
				   public G4CMPProcessUtils { 
public:
  G4CMPInterValleyScattering();
  virtual ~G4CMPInterValleyScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  G4CMPInterValleyScattering(G4CMPInterValleyScattering&);
  G4CMPInterValleyScattering& operator=(const G4CMPInterValleyScattering& right);

  G4VProcess* stepLimiter;
  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;

  G4double velLong;
  G4double l0Hole;
  G4double massHole;
  G4double ksound_Hole;
  G4double massFreeElectron;
  G4double hbar;
  G4double E_0_ED_203;
  G4double E_0_ED_201;
};

#endif
