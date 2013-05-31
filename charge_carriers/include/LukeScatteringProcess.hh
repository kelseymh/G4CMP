#ifndef LukeScatteringProcess_h
#define LukeScatteringProcess_h 1
#define PI 3.14153

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class LukeScatteringProcess : public G4VDiscreteProcess
{
public:


  G4double MakeTheta(G4double& k,G4double& ks);
  G4double MakePhi(G4double& k, G4double& ks, G4double& theta);

  //LukeScatteringProcess(const G4String& processName = "LukeScattering");
  LukeScatteringProcess();

  virtual ~LukeScatteringProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  LukeScatteringProcess(LukeScatteringProcess&);
  LukeScatteringProcess& operator=(const LukeScatteringProcess& right);

  G4VProcess* stepLimiter;

  G4double velLong;
  G4double l0Hole;
  G4double massHole;
  G4double ksound_Hole;
  G4double massFreeElectron;
  G4double hbar;

};



#endif
