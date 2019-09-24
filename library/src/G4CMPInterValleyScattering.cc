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
// 20180831  Change G4CMPInterValleyScattering to use Lin. and Quad. models
// 20190704  Add selection of rate model by name, and material specific
// 20190904  C. Stanford -- Add 50% momentum flip (see G4CMP-168)
// 20190906  Push selected rate model back to G4CMPTimeStepper for consistency

#include "G4CMPInterValleyScattering.hh"
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


// Construcor and destructor

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("G4CMPInterValleyScattering", fInterValleyScattering) {
  UseRateModel(G4CMPConfigManager::GetIVRateModel());
}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering() {;}


// Only electrons have physical valleys associated with them

G4bool 
G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsElectron(&aPD);
}


// Select different rate models by string (globally or by material)

void G4CMPInterValleyScattering::UseRateModel(G4String model) {
  if (model.empty()) {			// No argument, initialize w/global
    if (GetRateModel()) return;		// Do not change existing model
    model = G4CMPConfigManager::GetIVRateModel();
  }

  model.toLower();
  if (model == modelName) return;	// Requested model already in use

  // Select from valid names; fall back to Quadratic if invalid name specified
       if (model(0) == 'q') UseRateModel(new G4CMPIVRateQuadratic);
  else if (model(0) == 'l') UseRateModel(new G4CMPIVRateLinear);
  else if (model(0) == 'i') UseRateModel(new G4CMPInterValleyRate);
  else {
    G4cerr << GetProcessName() << " ERROR: Unrecognized rate model '"
	   << model << "'" << G4endl;
    if (!GetRateModel()) UseRateModel("Quadratic");
  }

  modelName = model;

  // Ensure that TimeStepper process is given new model
  PushModelToTimeStepper();
}


// Switch rate models if necessary based on material

G4double G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& track,
						     G4double prevStep,
						     G4ForceCondition* cond) {
  UseRateModel(theLattice->GetIVModel());	// Use current material's rate
  return G4CMPVProcess::GetMeanFreePath(track, prevStep, cond);
}


// Perform scattering action

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }
  
  // Don't do anything at a volume boundary
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
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
  
  // There's a 50% chance that the charge jumped into the antivalley rather
  // than the primary valley. If so, its momentum needs to be reversed to 
  // preserve symmetry.
  if (G4UniformRand()>0.5) p = -p;

  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);
  
  ClearNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}


// Ensure the same rate model is used here and in G4CMPTimeStepper

void G4CMPInterValleyScattering::PushModelToTimeStepper() {
  if (verboseLevel>1)
    G4cout << GetProcessName() << "::PushModelToTimeStepper" << G4endl;

  G4CMPTimeStepper* tsProc =
    dynamic_cast<G4CMPTimeStepper*>(G4CMP::FindProcess(GetCurrentTrack(),
						       "G4CMPTimeStepper"));

  if (tsProc) tsProc->UseIVRateModel(GetRateModel());
}
