/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140324  Restore Z-axis mass tensor
// 20140331  Add required process subtype code
// 20140418  Drop local valley transforms, use lattice functions instead
// 20140429  Recompute kinematics relative to new valley
// 20140908  Allow IV scatter to change momentum by conserving energy
// 20150109  Revert IV scattering to preserve momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160601  Must apply lattice rotation before valley.
// 20161004  Change valley selection function to avoid null choice
// 20161114  Use G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for particle identity checks
// 20170802  Remove MFP calculation; use scattering-rate model
// 20170809  Replace Edelweiss rate with physical (matrix element) model
// 20170821  Use configuration flag to choose Edelweiss vs. physical rate

#include "G4CMPInterValleyScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPInterValleyRate.hh"
#include "G4CMPIVRateEdelweiss.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include <math.h>

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("G4CMPInterValleyScattering", fInterValleyScattering) {
  if (G4CMPConfigManager::UseIVEdelweiss()) 
    UseRateModel(new G4CMPIVRateEdelweiss);
  else
    UseRateModel(new G4CMPInterValleyRate);
}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering() {;}


G4bool 
G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsElectron(&aPD);
}


G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  // Do nothing when step limit reached or leaving volume
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }

  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return &aParticleChange;
  }

  // Get track's energy in current valley
  G4ThreeVector p = GetLocalMomentum(aTrack);
  G4int valley = GetValleyIndex(aTrack);
  p = theLattice->MapPtoK_valley(valley, p); // p is actually k now

  // picking a new valley at random if IV-scattering process was triggered
  valley = ChangeValley(valley);
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);

  p = theLattice->MapK_valleyToP(valley, p); // p is p again
  RotateToGlobalDirection(p);

  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);

  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}
