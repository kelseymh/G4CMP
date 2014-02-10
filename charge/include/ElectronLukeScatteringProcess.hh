#ifndef ElectronLukeScatteringProcess_h
#define ElectronLukeScatteringProcess_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class ElectronLukeScatteringProcess : public G4VDiscreteProcess
{
public:


  G4double MakeTheta(G4double& k,G4double& ks);
  G4double MakePhi(G4double& k, G4double& ks, G4double& theta);

  //ElectronLukeScatteringProcess(const G4String& processName = "LukeScattering");
  ElectronLukeScatteringProcess(G4VProcess*);

  virtual ~ElectronLukeScatteringProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  ElectronLukeScatteringProcess(ElectronLukeScatteringProcess&);
  ElectronLukeScatteringProcess& operator=(const ElectronLukeScatteringProcess& right);

  G4AffineTransform NormalToValley;
  G4AffineTransform ValleyToNormal;
  //G4AffineTransform HVTransform;
  G4RotationMatrix mInv; // Inverse mass tensor
  G4ThreeVector T;       // HV Transformation matrix diagonal
  G4VProcess* stepLimiter;

  G4double velLong;      // Speed of sound of longitudinal phonons
  G4double l0;           // Characteristic scattering length
  G4double me;           // Mass of electron
  G4double mc;           // Effective mass of electron (units of electron)
  G4double ksound;       // Wave number for mc moving at velLong
};



#endif
