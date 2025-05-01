/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
 \***********************************************************************/

/// \file library/include/G4CMPIVSwitchModel.hh
/// \brief Simple but non-physical implementation of intervalley scattering
///        of electrons during charge transport.  Reassigns valley index
///        registered for track, moving the momentum direction to preserve
///	   the angle with respect to the valley axis (conserves energy).
//
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
// 20180831  Change G4CMPIVSwitchModel to use Lin. and Quad. models
// 20190704  Add selection of rate model by name, and material specific
// 20190904  C. Stanford -- Add 50% momentum flip (see G4CMP-168)
// 20190906  Push selected rate model back to G4CMPTimeStepper for consistency
// 20231122  Remove 50% momentum flip (see G4CMP-375)
// 20240823  Allow ConfigManager IVRateModel setting to override config.txt
// 20250430  Move PostStepDoIt() implementation from G4CMPInterValleyScattering

#include "G4CMPIVSwitchModel.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPInterValleyRate.hh"
#include "G4CMPIVRateQuadratic.hh"
#include "G4CMPIVRateLinear.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include <math.h>


// Move electron from one valley to another without changing energy

G4VParticleChange* G4CMPIVSwitchModel::PostStepDoIt(const G4Track& aTrack, 
						    const G4Step& aStep) {
  InitializeParticleChange(GetValleyIndex(aTrack), aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  if (verboseLevel > 1) {
    G4cout << modelName << "::PostStepDoIt: Step limited by "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }
  
  // Don't do anything at a volume boundary
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return &aParticleChange;
  }
  
  // Get track's energy in current valley
  G4ThreeVector p = GetLocalMomentum(aTrack);
  G4int valley = GetValleyIndex(aTrack);
  p = theLattice->MapPtoK(valley, p); // p is actually k now
  p = theLattice->RotateToValley(valley, p);
  
  // picking a new valley at random if IV-scattering process was triggered
  valley = ChangeValley(valley);
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);

  p = theLattice->RotateFromValley(valley, p);
  p = theLattice->MapKtoP(valley, p); // p is p again
  RotateToGlobalDirection(p);
  
  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);
  
  return &aParticleChange;
}
