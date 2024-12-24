/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBogoliubovQPRandomWalkBoundary.hh
/// \brief Definition of the  G4CMPBogoliubovQPRandomWalkBoundary class

#ifndef G4CMPBogoliubovQPRandomWalkBoundary_h
#define G4CMPBogoliubovQPRandomWalkBoundary_h 1

#include "G4VBogoliubovQPProcess.hh"
#include "G4CMPSCUtils.hh"
#include "G4CMPBoundaryUtils.hh"

class G4CMPProcessUtils;

class G4CMPBogoliubovQPRandomWalkBoundary : public G4VBogoliubovQPProcess,
					    public G4CMPBoundaryUtils {
public:
  G4CMPBogoliubovQPRandomWalkBoundary(const G4String& processName="G4CMPBogoliubovQPRandomWalkBoundary");
  
  virtual ~G4CMPBogoliubovQPRandomWalkBoundary();

  //  // Configure for current track including AnharmonicDecay utility
  //  virtual void LoadDataForTrack(const G4Track* track);
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
							G4double previousStepSize,
							G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);
                    
                       
  G4bool CheckQPVolumes(const G4Step& aStep);
  G4bool IsValidQPVolume(G4VPhysicalVolume* volume,
			 G4double qpEKin );

protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double prevStepLength,
                                   G4ForceCondition* condition);

  // Apply QP-specific conditions, after calling through to base
  virtual G4bool ReflectTrack(const G4Track& aTrack, const G4Step& aStep) const;
  virtual void DoAbsorption(const G4Track& aTrack, const G4Step& aStep, G4ParticleChange& aParticleChange);
  virtual void DoReflection(const G4Track& aTrack, const G4Step& aStep, G4ParticleChange& aParticleChange);
  virtual void DoTransmission(const G4Track& aTrack, const G4Step& aStep, G4ParticleChange& aParticleChange);

  G4ThreeVector GetLambertianVector(const G4ThreeVector& surfNorm) const;
  
  //Boolean to indicate whether the pre/post-step volumes have valid material properties for QP transport
  G4bool preQPVolume;
  G4bool postQPVolume;
   
  G4String procName;
  G4CMPProcessUtils* procUtils;
   
private:

  // hide assignment operator as private
  G4CMPBogoliubovQPRandomWalkBoundary(G4CMPBogoliubovQPRandomWalkBoundary&);
  G4CMPBogoliubovQPRandomWalkBoundary(G4CMPBogoliubovQPRandomWalkBoundary&&);
  G4CMPBogoliubovQPRandomWalkBoundary& operator=(const G4CMPBogoliubovQPRandomWalkBoundary&);
  G4CMPBogoliubovQPRandomWalkBoundary& operator=(const G4CMPBogoliubovQPRandomWalkBoundary&&);
};

#endif	/* G4CMPBogoliubovQPRandomWalkBoundary_h */
