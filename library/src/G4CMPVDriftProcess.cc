/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVDriftProcess.cc
/// \brief Implementation of the G4CMPVDriftProcess base class
//
// $Id$
//
// 20140325  Move time-step calculation here from TimeStepper and LukeScat
// 20140331  Add required process subtype code
// 20141216  Set true electron velocity (and force it) in SetNewKinematics()
// 20150109  Use G4CMP_SET_ELECTRON_MASS to enable dynamical mass, velocity
// 20150111  Add functionality to enforce minimum step length
// 20150112  Handle holes as well as electrons with FillParticleChange()
// 20160624  Use GetTrackInfo() accessor
// 20160829  Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
// 20161114  Use new G4CMPDriftTrackInfo
// 20170601  Inherit from new G4CMPVProcess, which provides G4CMPProcessUtils
// 20170620  Follow interface changes in G4CMPProcessUtils
// 20201231  FillParticleChange() should also reset valley index if requested

#include "G4CMPVDriftProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessType.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"


// Constructor and destructor
// NOTE:  Initial values are arbitrary and non-physical
G4CMPVDriftProcess::G4CMPVDriftProcess(const G4String& processName,
                                       G4CMPProcessSubType stype)
  : G4CMPVProcess(processName, stype) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPVDriftProcess::~G4CMPVDriftProcess() {;}


// Only applies to the known charge carriers

G4bool G4CMPVDriftProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsChargeCarrier(aPD);
}


// Overload base version to set a minimum step size, avoiding "stuck" tracks

G4double
G4CMPVDriftProcess::PostStepGetPhysicalInteractionLength(
                      const G4Track& track,
                      G4double previousStepSize,
                      G4ForceCondition* condition) {
  G4double trueLength =
    G4VDiscreteProcess::PostStepGetPhysicalInteractionLength(track,
                                                             previousStepSize,
                                                             condition);

  G4double minLength = G4CMPConfigManager::GetMinStepScale();
  minLength *= (IsElectron() ? theLattice->GetElectronScatter()
		: theLattice->GetHoleScatter());

  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepGPIL: minLength " << minLength
	   << " trueLength " << trueLength << G4endl;
  }
  
  return minLength<trueLength ? trueLength : minLength;
}


// Fill ParticleChange energy and mass for electron charge carrier momentum

void 
G4CMPVDriftProcess::FillParticleChange(G4int ivalley, const G4ThreeVector& p) {
  G4double mass = GetCurrentTrack()->GetDynamicParticle()->GetMass();

  G4ThreeVector v;
  if (IsElectron()) {
    v = GetGlobalDirection(theLattice->MapPtoV_el(ivalley, GetLocalDirection(p)));
  } else if (IsHole()) {
    v = p*c_light/mass;		// p and mass in MeV, not MeV/c, MeV/c^2
  } else {
    G4Exception("G4CMPVDriftProcess::FillParticleChange", "DriftProcess001",
    EventMustBeAborted, "Unknown charge carrier");
  }

  // Non-relativistic, but "mass" is mc^2
  G4double energy = 0.5*mass*v.mag2()/c_squared;

  FillParticleChange(ivalley, energy, v);
}

// Fill ParticleChange mass for electron charge carrier with given energy

void 
G4CMPVDriftProcess::FillParticleChange(G4int ivalley, G4double Ekin,
             const G4ThreeVector& v) {
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(GetCurrentTrack())->SetValleyIndex(ivalley);

  aParticleChange.ProposeMomentumDirection(v.unit());
  aParticleChange.ProposeEnergy(Ekin);
}
