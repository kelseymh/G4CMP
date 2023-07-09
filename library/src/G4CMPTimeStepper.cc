/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils, also
//	     use proper TransformAxis() on vectors, *not* TransformPoint()
//	     Add wrapper function to compute individual time steps
// 20140314  Fetch propagation parameters from lattice, instead of hardwired
// 20140324  Migrate to use of volume-local field, do coordinate transforms
// 20140325  Move most of time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140418  Remove explicit valley transforms, use lattice function
// 20140429  Adjust "effective mass" and energy based on end-of-step momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange, drop
//	     redundant IsApplicable()
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160624  Use GetTrackInfo() accessor
// 20161114  Use new G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for track identity functions
// 20170806  Swap GPIL and MFP functions to work with G4CMPVProcess base
// 20170905  Cache Luke and IV rate models in local LoadDataFromTrack()
// 20170908  Remove "/10." rescaling of field when computing steps
// 20170919  Use rate threshold interface to define alternate step lengths
// 20180831  Fix compiler warning with PostStepDoIt() arguments
// 20190906  Provide functions to externally set rate models, move process
//		lookup functionality to G4CMP(Track)Utils.
// 20200331  C. Stanford (G4CMP-195): Added charge trapping
// 20200331  G4CMP-196: Added impact ionization mean free path
// 20200426  G4CMP-196: Use static function in TrapIonization for MFPs
// 20200504  M. Kelsey (G4CMP-195):  Get trapping MFPs from process
// 20200520  "First report" flag must be thread-local.
// 20200804  Move field access to G4CMPFieldUtils
// 20210923  Ensure that rate calculations are initialized for track
// 20230527  Drop competitive MFP calculation; use this process to enforce
//	       a maximum allowed step length, to support mass recalculation.

#include "G4CMPTimeStepper.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {;}


// Return configured maximum step length, or DBL_MAX if not set

G4double G4CMPTimeStepper::GetMeanFreePath(const G4Track& aTrack, G4double,
					   G4ForceCondition* cond) {
  *cond = NotForced;

  G4double maxStep = G4CMPConfigManager::GetMaximumStep();

  G4ThreadLocal static G4bool first=true;
  if (verboseLevel>1 && first) {
    G4cout << "G4CMPTimeStepper::GetMFP using maxStep " << maxStep << G4endl;
    first = false;
  }
  
  return (maxStep>0. ? maxStep : DBL_MAX);
}


// At end of step, recompute kinematics; important for electrons

G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  // Adjust dynamical mass for electrons using end-of-step momentum direction
  G4ThreeVector plocal = GetLocalMomentum(aTrack);
  G4double meff = IsHole() ? theLattice->GetHoleMass()
    : theLattice->GetElectronEffectiveMass(GetValleyIndex(aTrack), plocal);

  if (IsElectron()) aParticleChange.ProposeMass(meff*c_squared);

  // Report basic kinematics
  if (verboseLevel) {
    G4cout << GetProcessName() << (IsElectron()?" elec":" hole")
	   << " Ekin " << GetKineticEnergy(aTrack)/eV << " eV,"
	   << " m " << meff*c_squared/electron_mass_c2 << " m_e," << G4endl
	   << " p " << GetGlobalMomentum(aTrack)/eV << " "
	   << GetGlobalMomentum(aTrack).mag()/eV << " eV"
	   << G4endl;
  }
  
  // Report electric field info (not valid if LukeScattering enabled)
  if (verboseLevel>1) {
    const G4StepPoint* pre = aStep.GetPreStepPoint();
    const G4StepPoint* post = aStep.GetPostStepPoint();
    G4double E0 = pre->GetKineticEnergy();
    G4double Ef = post->GetKineticEnergy();
    G4ThreeVector p0 = pre->GetMomentum();
    G4ThreeVector pf = post->GetMomentum();
    G4ThreeVector pos0 = pre->GetPosition();
    G4ThreeVector posf = post->GetPosition();

    G4ThreeVector field = G4CMP::GetFieldAtPosition(aTrack.GetTouchable(),pos0);
    G4ThreeVector deltaP = pf-p0;
    G4double deltaE = Ef-E0;
    G4double deltaV = field.dot(posf-pos0);	// Voltage drop

    G4cout << " pre-step @ " << pos0 << G4endl << "   E " << E0/eV << " eV,"
	   << " p0 " << p0/eV << " " << p0.mag()/eV << " eV" << G4endl
	   << " post-step @ " << posf << G4endl << "   E " << Ef/eV << " eV,"
	   << " pf " << pf/eV << " " << pf.mag()/eV << " eV" << G4endl
	   << " E-field " << field/(volt/m) << " " << field.mag()/(volt/m)
	   << " V/m" << G4endl
	   << " dPos " << posf-pos0 << " " << (posf-pos0).mag()/mm << " mm,"
	   << G4endl
	   << " dV = field.dPos " << deltaV/volt << " V,"
	   << " dE " << deltaE/eV << " eV" << G4endl
	   << " dP " << deltaP/eV << " " << deltaP.mag()/eV << " eV"
	   << G4endl;

    G4double EvsV = deltaE - (IsElectron()?-1:1)*deltaV;
    if (fabs(EvsV) > 1e-12) {
      G4cout << "*** Energy-voltage mismatch: |dE-qdV| " << EvsV/eV
	     << " > 1e-6 eV" << G4endl;
    }
  }

  ClearNumberOfInteractionLengthLeft();		// All processes must do this!
  return &aParticleChange;
}
