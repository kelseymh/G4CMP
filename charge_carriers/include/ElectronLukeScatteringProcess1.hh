#ifndef ElectronLukeScatteringProcess1_h
#define ElectronLukeScatteringProcess1_h 1
#define PI 3.14153

#include "globals.hh"
#include "G4VProcess.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class ElectronLukeScatteringProcess1 : public G4VProcess
{
public:


  G4double MakeTheta(G4double& k,G4double& ks);
  G4double MakePhi(G4double& k, G4double& ks, G4double& theta);

  //ElectronLukeScatteringProcess(const G4String& processName = "LukeScattering");
  ElectronLukeScatteringProcess1();

  virtual ~ElectronLukeScatteringProcess1();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:

 virtual G4double PostStepGetPhysicalInteractionLength( 
						const G4Track& aTrack,
						G4double prevStepSize,
						G4ForceCondition* cond);


  //no-op in atRest GPIL
  virtual G4double AtRestGetPhysicalInteractionLength(
						      const G4Track&,
						      G4ForceCondition*
						      ){ return -1.0;};

  //no-op in AtRestDoIt
  virtual G4VParticleChange* AtRestDoIt(
					const G4Track&,
					const G4Step&
					){ return 0;};

  //no-op in alongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(
						      const G4Track&,
						      G4double,
						      G4double,
						      G4double&,
						      G4GPILSelection*
						      ){ return -1.0;};

  //no-op in alongStepDoIt
  virtual G4VParticleChange* AlongStepDoIt(
					const G4Track&,
					const G4Step&
					){ return 0;};
private:
  //hide assignment operator as private
  ElectronLukeScatteringProcess1(ElectronLukeScatteringProcess1&);
  ElectronLukeScatteringProcess1& operator=(const 
ElectronLukeScatteringProcess1& right);

  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;
  //G4AffineTransform HVTransform;
  //G4RotationMatrix HVMatrix;
  G4VProcess* stepLimiter;

  G4double velLong;
  G4double l0Hole;
  G4double massHole;
  G4double ksound_Hole;
  G4double massFreeElectron;
  G4double hbar;





};



#endif
