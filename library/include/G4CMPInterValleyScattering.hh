// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice

#ifndef G4CMPInterValleyScattering_h
#define G4CMPInterValleyScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class G4CMPInterValleyScattering : public G4CMPVDriftProcess { 
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
};

#endif
