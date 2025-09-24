/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPQPDiffusion.hh
/// \brief Definition of the G4CMPQPDiffusion class
///
/// This is the class definition for the process that does QP diffusion
/// using the walk on spheres technique. For more detailed info, see
/// the class implementation.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPQPDiffusion_h
#define G4CMPQPDiffusion_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPSCUtils.hh"
#include "G4CMPBoundaryUtils.hh"
#include "G4CMPParticleChangeForQPDiffusion.hh"
#include "G4CMPProcessSubType.hh"
#include "G4ParticleChange.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include <fstream> 

class G4ParticleDefinition;
class G4SafetyHelper;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class G4CMPQPDiffusion : public G4VContinuousDiscreteProcess,
			 public G4CMPProcessUtils, public G4CMPSCUtils {
public:
  G4CMPQPDiffusion(const G4String& name="qpDiffusion",
		   G4CMPProcessSubType stype=fQPDiffusion);

  virtual ~G4CMPQPDiffusion();

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

  G4ThreeVector
  FindDirectionToNearbyBoundary(const G4Track& track,
				const G4ThreeVector& trackPosition,
				const G4double the2DSafety,
				G4bool & needToRepeatCalculation,
				G4bool useSweepForDaughterSafety=false);

  G4bool IsApplicable(const G4ParticleDefinition& aPD);
  
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
  G4CMPQPDiffusion(G4CMPQPDiffusion &) = delete;
  G4CMPQPDiffusion &
    operator=(const G4CMPQPDiffusion &right) = delete;
    
  // ======== Parameters of the class fixed at initialisation =======
  //Safety helper to query G4Navigator and check distance to geometric
  //boundaries
  G4SafetyHelper*             fSafetyHelper;
  
  G4double SampleTimeStepFromFirstPassageDistribution(G4double the2DSafety);
  G4double SampleDimensionlessTimeStepUsingAcceptanceRejectionTechnique();
  
protected:
  virtual G4bool UpdateMeanFreePathForLatticeChangeover(const G4Track& aTrack);
  virtual void UpdateSCAfterLatticeChange();
  
  //Custom particle change class where the UpdateAlongStep method handles
  //non-physical changes to particle from random walk process
  G4CMPParticleChangeForQPDiffusion fParticleChange;
    
private:
  G4double      fTimeStep;               //Time increment for the step
  G4double      fPathLength;             //Path length returned by the AlongStepGPIL (starts diffusion-unfolded, and then folds in diffusion)
  G4double      fPreDiffusionPathLength; //Initial, diffusion-unfolded path length, for persistency
  G4double      fDiffConst;              //Energy dependent diffusion constant
  G4double      f2DSafety;               //The 2D safety computed for this step
  G4double      fTimeStepToBoundary;     //Timestep computed for step's 2Dsafety
  G4ThreeVector fOldPosition;            //Position at the beginning of the step
  G4ThreeVector fNewPosition;            //Proposed position after diffusion
  G4ThreeVector fNewDirection;           //Direction of proposed diffusion step 
  G4bool        fPositionChanged=false;
  G4bool        isActive=false;
  G4bool        fTrackOnBoundary=false;
  G4bool        fVerySmallStep=false;
  G4bool        fVerySmallStepInsideSoftFloor=false;
  G4double      fBoundaryFudgeFactor;
  G4double      fHardFloorBoundaryScale;  
  G4double      fSoftFloorBoundaryScale; //Boundary eps used in walkonspheres
  G4double      fPicometerScale;
  
  //The last N boundary scatters for this
  std::vector<std::pair<G4ThreeVector,G4ThreeVector> > fBoundaryHistory;
  
  G4int fBoundaryHistoryTrackID;           //Which QP to produce bdry history of
  G4int fMaxBoundaryHistoryEntries;        //Last N QP boundary scatters
  G4bool fQPIsStuck;                       //Is QP stuck in corner?
  G4ThreeVector fStuckNorm1;               //For QPs stuck in corner
  G4ThreeVector fStuckNorm2;               //For QPs stuck in corner
  G4ThreeVector fOutgoingSurfaceTangent1;  //For QPs stuck in corner
  G4ThreeVector fOutgoingSurfaceTangent2;  //For QPs stuck in corner
  G4double fDotProductDefiningUniqueNorms; //Self explanatory
  G4double fStuckInCornerThreshold;        //How spatially tight until QP stuck?
  G4double fNeedSweptSafetyInGetMFP;       //Flag to redo safeties w/more care
  G4double fPreemptivelyKillTrack;         //Kill track
  
  void UpdateBoundaryHistory(G4int trackID, G4ThreeVector preStepPos,
			     G4ThreeVector preStepNorm);
  std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector>
  CheckForStuckQPs();
  
  std::tuple<G4bool,G4ThreeVector,G4ThreeVector,G4ThreeVector,G4ThreeVector>
  CheckForStuckQPsInCorner();
  
  std::pair<G4ThreeVector,G4ThreeVector>
  FindSurfaceTangentsForStuckQPEjection(G4ThreeVector norm1,
					G4ThreeVector pos1,
					G4ThreeVector norm2,
					G4ThreeVector pos2,
					G4ThreeVector & cornerLocation);
  
  void PostCheckBulkTreatment(G4double stepTransportOnlyDeltaT);

  G4bool CheckForPhantomBoundaryCrossings(G4ThreeVector trackPosition,
					  G4double the2DSafety,
					  G4double originalOption1Safety,
					  G4double originalOption2Safety,
					  G4ThreeVector outputDir);
  
  G4double HandleVerySmallSteps(G4double thisMFP, G4double the2DSafety,
				G4double velocity);
  G4double ComputePathLengthInGoldilocksZone(); 
};

// ======== Run time inline methods ================

inline G4double G4CMPQPDiffusion::TimeStep() const
{
  return fTimeStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CMPQPDiffusion::SetTimeStep(G4double val)
{
  fTimeStep = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
