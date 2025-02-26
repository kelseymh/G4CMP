/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPBogoliubovQPRandomWalkTransport.hh
/// \brief Definition of the G4CMPBogoliubovQPRandomWalkTransport class

#ifndef G4CMPBogoliubovQPRandomWalkTransport_h
#define G4CMPBogoliubovQPRandomWalkTransport_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4CMPParticleChangeForBogoliubovQPRandomWalk.hh"
#include "G4ParticleChange.hh"
#include "G4CMPProcessUtils.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4CMPProcessSubType.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4CMPSCUtils.hh"
#include "G4CMPBoundaryUtils.hh"

class G4ParticleDefinition;
class G4SafetyHelper;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4CMPBogoliubovQPRandomWalkTransport : public G4VContinuousDiscreteProcess, public G4CMPProcessUtils, public G4CMPBoundaryUtils, public G4CMPSCUtils
{
public:

    G4CMPBogoliubovQPRandomWalkTransport(const G4String& name = "G4CMPBogoliubovQPRandomWalkTransport", G4CMPProcessSubType stype = fBogoliubovQPRandomWalkTransport);

  virtual ~G4CMPBogoliubovQPRandomWalkTransport();

public:

  //------------------------------------------------------------------------
  // Generic methods common to all ContinuousDiscrete processes
  //------------------------------------------------------------------------

  // This is called in the beginning of tracking for a new track
  void StartTracking(G4Track*) override;

  // The function overloads the corresponding function of the base
  // class.It limits the step near to boundaries only
  // and invokes the method GetMscContinuousStepLimit at every step.
  G4double AlongStepGetPhysicalInteractionLength(
                                        const G4Track&,
                                        G4double  previousStepSize,
                                        G4double  currentMinimalStep,
                                        G4double& currentSafety,
                                        G4GPILSelection* selection) override;

  // The function overloads the corresponding function of the base
  // class.
  G4double PostStepGetPhysicalInteractionLength(
                                      const G4Track& track,
                                      G4double  previousStepSize,
                                      G4ForceCondition* condition) override;

  // Along step actions
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  // Post step actions
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  G4ThreeVector FindDirectionToNearbyBoundary(const G4Track& track, const G4double the2DSafety );
  
public:
 
  // Update the default fTimeStep of the random walker if currentMinimalStep is longer
  inline G4double TimeStep() const;
  inline void SetTimeStep(G4double val);

  //------------------------------------------------------------------------
  // Run time methods
  //------------------------------------------------------------------------

protected:

  // This method is not used for tracking, it returns mean free path value
  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  // This method is not used for tracking, it returns step limit
  G4double GetContinuousStepLimit(const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimalStep,
                                  G4double& currentSafety) override ;

  G4double ContinuousStepLimit(const G4Track& track,
                               G4double previousStepSize,
                               G4double currentMinimalStep,
                               G4double& currentSafety);

private:

  // hide  assignment operator
  G4CMPBogoliubovQPRandomWalkTransport(G4CMPBogoliubovQPRandomWalkTransport &) = delete;
  G4CMPBogoliubovQPRandomWalkTransport &
    operator=(const G4CMPBogoliubovQPRandomWalkTransport &right) = delete;
    
  // ======== Parameters of the class fixed at initialisation =======
  //Safety helper to query G4Navigator and check distance to geometric boundaries
  G4SafetyHelper*             fSafetyHelper;

  // ======== Cached values - may be state dependent ================

protected:

  virtual G4bool UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack);
  virtual void UpdateSCAfterLatticeChange();
  
  //Custom particle change class where the UpdateAlongStep method handles
  //non-physical changes to particle from random walk process
  G4CMPParticleChangeForBogoliubovQPRandomWalk fParticleChange;
    
private:
    
  G4double                    fTimeStep;               //Time increment for the step
  G4double                    fPathLength;             //Path length returned by the AlongStepGPIL (starts diffusion-unfolded, and then folds in diffusion)
  G4double                    fPreDiffusionPathLength; //Initial, diffusion-unfolded path length, for persistency
  G4double                    fDiffConst;              //Energy dependent diffusion constant

  G4double                    f2DSafety;               //The 2D safety computed for this step
  G4double                    fTimeStepToBoundary;     //The time step computed given the 2D safety computed for this step
  
  G4ThreeVector               fOldPosition;            //Position at the beginning of the step
  G4ThreeVector               fNewPosition;            //Proposed position after diffusion
  G4ThreeVector               fNewDirection;           //Direction of proposed effective diffusion step
    
  G4bool                      fPositionChanged= false;
  G4bool                      isActive= false;
  G4bool                      fTrackOnBoundary= false;
  G4bool                      fVerySmallStep= false;
  G4double                    fBoundaryFudgeFactor;
};

// ======== Run time inline methods ================

inline G4double G4CMPBogoliubovQPRandomWalkTransport::TimeStep() const
{
  return fTimeStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CMPBogoliubovQPRandomWalkTransport::SetTimeStep(G4double val)
{
  fTimeStep = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
