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
// 20230527  G4CMP-295: Adjust base class interaction length parameters if
//	        minimum path length used.

#include "G4CMPVDriftProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPVScatteringRate.hh"
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
#include <algorithm>
#include <initializer_list>

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
  // Get desired step length computed from MFP and number of ILs
  G4double trueLength =
    G4VDiscreteProcess::PostStepGetPhysicalInteractionLength(track,
                                                             previousStepSize,
                                                             condition);

  // Check for upcoming energy threshold in rate
  G4double ekin = GetKineticEnergy(track);
  G4CMPVScatteringRate* processRate = GetRateModel();
  G4double energyStep = (processRate ? EnergyStep(processRate->Threshold(ekin))
			 : -1.);

  // Get configured minimum step (zero if not configured)
  G4double minLength = G4CMPConfigManager::GetMinimumStep();

  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepGPIL:"
	   << " trueLength " << trueLength/mm << " mm"
	   << " minLength " << minLength/mm << " mm"
	   << " energyStep " << energyStep/mm << " mm" << G4endl;
  }

  // If desired step is shorter than cutoff, force process now
  if (minLength > 0. && trueLength < minLength) {
    theNumberOfInteractionLengthLeft = 0.;	// Assert that we've reached end
    return 0.;
  }

  // If threshold happens before desired step, override IL, #IL, step length
  if (energyStep > 0. && energyStep < trueLength) {
    theNumberOfInteractionLengthLeft = 1.;
    currentInteractionLength = energyStep;
    return energyStep;
  }

  return trueLength;		// Nothing special, use Geant4 step
}


// Fill ParticleChange energy and mass for electron charge carrier momentum

void 
G4CMPVDriftProcess::FillParticleChange(G4int ivalley, const G4ThreeVector& p) {
  // Compute kinetic energy from momentum for electrons or holes
  G4double energy =
    (IsElectron() ? theLattice->MapPtoEkin(ivalley, GetLocalDirection(p))
     : p.mag2()/(2.*GetCurrentTrack()->GetDynamicParticle()->GetMass()) );

  FillParticleChange(ivalley, energy, p);
}

// Fill ParticleChange mass for electron charge carrier with given energy

void G4CMPVDriftProcess::FillParticleChange(G4int ivalley, G4double Ekin,
					    const G4ThreeVector& v) {
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(GetCurrentTrack())->SetValleyIndex(ivalley);

  aParticleChange.ProposeMomentumDirection(v.unit());
  aParticleChange.ProposeEnergy(Ekin);

  if (IsElectron()) {		// Geant4 wants mc^2, not plain mass
    G4double meff = theLattice->GetElectronEffectiveMass(ivalley,
							 GetLocalDirection(v));
    aParticleChange.ProposeMass(meff*c_squared);
  }
}

G4double G4CMPVDriftProcess::EnergyStep(G4double Efinal) const {
  const G4Track* trk = GetCurrentTrack();

  // NOTE: Assumes particle travels along local field line
  G4double EField = G4CMP::GetFieldAtPosition(*trk).mag();
  if (EField <= 0.) return -1.;			// No field, no acceleration

  G4double Ekin = GetKineticEnergy(trk);
  if (Ekin > Efinal) return -1.;		// Already over threshold

  if (verboseLevel>1) {
    G4cout << "G4CMPVDriftProcess::EnergyStep from " << Ekin/eV
	   << " to " << Efinal/eV << " eV" << G4endl;
  }

  // Add 20% rescaling to account for electron oblique propagation
  return 1.2*(Efinal-Ekin)/EField;
}
