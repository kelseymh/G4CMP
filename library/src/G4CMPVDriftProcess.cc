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
  : G4CMPVProcess(processName, stype), velLong(330*m/s) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPVDriftProcess::~G4CMPVDriftProcess() {;}


// Only applies to the known charge carriers

G4bool G4CMPVDriftProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsChargeCarrier(&aPD);
}


// Get additional parameters from lattice for carriers

void G4CMPVDriftProcess::LoadDataForTrack(const G4Track* track) {
  G4CMPProcessUtils::LoadDataForTrack(track);
  velLong = theLattice->GetSoundSpeed();
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

  const G4double scale = G4CMPConfigManager::GetMinStepScale();

  if (scale > 0.) {
    G4double l0 = 0.;
    if (IsElectron(&track)) {
      l0 = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(track)->Lattice()->GetElectronScatter();
    } else if (IsHole(&track)) {
      l0 = G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(track)->Lattice()->GetHoleScatter();
    } else {
      G4Exception("G4CMPVDriftProcess::FillParticleChange", "DriftProcess002",
                  EventMustBeAborted, "Unknown charge carrier");
    }

    G4double minLength = scale * l0;

    if (verboseLevel > 1) {
      G4cout << GetProcessName() << "::PostStepGPIL: minLength " << minLength
	     << " trueLength " << trueLength << G4endl;
    }

    return minLength<trueLength ? trueLength : minLength;
  }

  return trueLength;
}


// Compute characteristic time step for charge carrier
// Parameters are "Mach number" (ratio with sound speed) and scattering length

G4double 
G4CMPVDriftProcess::ChargeCarrierTimeStep(G4double mach, G4double l0) const {
  if (mach < 1.) return 3*l0/velLong;	// Sanity check if below sound speed

  return ( 3*l0 * mach / (velLong * (mach-1)*(mach-1)*(mach-1)) );
}


// Fill ParticleChange energy and mass for electron charge carrier momentum

void 
G4CMPVDriftProcess::FillParticleChange(G4int ivalley, const G4ThreeVector& p) {
  G4double mass = GetCurrentTrack()->GetDynamicParticle()->GetMass();

  G4ThreeVector v;
  if (GetCurrentParticle() == G4CMPDriftElectron::Definition()) {
    G4ThreeVector p_local = G4CMP::GetLocalDirection(GetCurrentVolume(), p);
    v = G4CMP::GetGlobalDirection(GetCurrentVolume(),
                                  theLattice->MapPtoV_el(ivalley, p_local));
  } else if (GetCurrentParticle() == G4CMPDriftHole::Definition()) {
    v = p*c_light/mass;
  } else {
    G4Exception("G4CMPVDriftProcess::FillParticleChange", "DriftProcess001",
    EventMustBeAborted, "Unknown charge carrier");
  }

  // Non-relativistic, but "mass" is mc^2
  G4double energy = 0.5*mass*v.mag2()/c_squared;

  FillParticleChange(ivalley, energy, v);
}

// Fill ParticleChange mass for electron charge carrier with given energy

// FIXME: Remove this trivial method, or at least remove the first parameter
// from its signature.
void 
G4CMPVDriftProcess::FillParticleChange(G4int /*ivalley*/, G4double Ekin,
             const G4ThreeVector& v) {
  aParticleChange.ProposeMomentumDirection(v.unit());
  aParticleChange.ProposeEnergy(Ekin);
}
