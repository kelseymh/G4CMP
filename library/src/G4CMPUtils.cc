/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPUtils.cc
/// \brief Namespace for general purpose or static utilities supporting
///	G4CMP.  Functions not dependent on current track/step info will
///     be moved here from G4CMPProcessUtils.
//
// $Id$
//
// 20170602  Provide call-by-reference versions of track identity functions
// 20170802  Provide scale factor argument to ChooseWeight functions
// 20170928  Replace "polarization" with "mode"
// 20190906  M. Kelsey -- Add function to look up process for track
// 20220816  M. Kelsey -- Move RandomIndex here for more general use
// 20220921  G4CMP-319 -- Add utilities for thermal (Maxwellian) distributions
// 20241223  G4CMP-419 -- Add utility to create per-thread debugging file
// 20250130  G4CMP-453 -- Apply coordinate rotations in PhononVelocityIsInward
// 20250422  G4CMP-468 -- Add displaced point test to PhononVelocityIsInward.
// 20250423  G4CMP-468 -- Add function to get diffuse reflection vector.
// 20250510  G4CMP-483 -- Ensure backwards compatibility for vector utilities.

#include "G4CMPUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4EventManager.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include <string>


// Select phonon mode using density of states in material

G4int G4CMP::ChoosePhononPolarization(const G4LatticePhysical* lattice) {
  return ChoosePhononPolarization(lattice->GetLDOS(),
                                  lattice->GetSTDOS(),
                                  lattice->GetFTDOS());
}

G4int G4CMP::ChoosePhononPolarization(G4double Ldos,
                                      G4double STdos, 
                                      G4double FTdos) {
  G4double norm = Ldos + STdos + FTdos;
  G4double cProbST = STdos/norm;
  G4double cProbFT = FTdos/norm + cProbST;

  // NOTE:  Order of selection done to match previous random sequences
  G4double modeMixer = G4UniformRand();
  if (modeMixer<cProbST) return G4PhononPolarization::TransSlow;
  if (modeMixer<cProbFT) return G4PhononPolarization::TransFast;
  return G4PhononPolarization::Long;
}

G4int G4CMP::ChooseValley(const G4LatticePhysical* lattice) {
  return static_cast<G4int>(G4UniformRand()*lattice->NumberOfValleys());
}


// Identify G4CMP particle categories

G4bool G4CMP::IsPhonon(const G4Track& track) {
  return IsPhonon(track.GetParticleDefinition());
}

G4bool G4CMP::IsPhonon(const G4Track* track) {
  return (track!=0 && IsPhonon(*track));
}

G4bool G4CMP::IsPhonon(const G4ParticleDefinition& pd) {
  return IsPhonon(&pd);
}

G4bool G4CMP::IsPhonon(const G4ParticleDefinition* pd) {
  return (G4PhononPolarization::Get(pd) != G4PhononPolarization::UNKNOWN);
}

G4bool G4CMP::IsElectron(const G4Track& track) {
  return IsElectron(track.GetParticleDefinition());
}

G4bool G4CMP::IsElectron(const G4Track* track) {
  return (track!=0 && IsElectron(*track));
}

G4bool G4CMP::IsElectron(const G4ParticleDefinition& pd) {
  return IsElectron(&pd);
}

G4bool G4CMP::IsElectron(const G4ParticleDefinition* pd) {
  return (pd == G4CMPDriftElectron::Definition());
}

G4bool G4CMP::IsHole(const G4Track& track) {
  return IsHole(track.GetParticleDefinition());
}

G4bool G4CMP::IsHole(const G4Track* track) {
  return (track!=0 && IsHole(*track));
}

G4bool G4CMP::IsHole(const G4ParticleDefinition& pd) {
  return IsHole(&pd);
}

G4bool G4CMP::IsHole(const G4ParticleDefinition* pd) {
  return (pd == G4CMPDriftHole::Definition());
}

G4bool G4CMP::IsChargeCarrier(const G4Track& track) {
  return (IsElectron(track) || IsHole(track));
}

G4bool G4CMP::IsChargeCarrier(const G4Track* track) {
  return (IsElectron(track) || IsHole(track));
}

G4bool G4CMP::IsChargeCarrier(const G4ParticleDefinition& pd) {
  return (IsElectron(pd) || IsHole(pd));
}

G4bool G4CMP::IsChargeCarrier(const G4ParticleDefinition* pd) {
  return (IsElectron(pd) || IsHole(pd));
}


// Generate weighting factor for phonons, charge carriers
// NOTE:  biasScale < 0. means to use "primary generator" scaling
// NOTE:  If zero is returned, track should NOT be created!

G4double G4CMP::ChooseWeight(const G4ParticleDefinition* pd,
			     G4double biasScale) {
  return (IsChargeCarrier(pd) ? ChooseChargeWeight(biasScale)
	  : IsPhonon(pd) ? ChoosePhononWeight(biasScale) : 1.);
}

G4double G4CMP::ChoosePhononWeight(G4double prob) {
  if (prob < 0.) prob = G4CMPConfigManager::GetGenPhonons();

  // If prob=0., random throw always fails, never divides by zero
  return ((prob==1.) ? 1. : (G4UniformRand()<prob) ? 1./prob : 0.);
}

G4double G4CMP::ChooseChargeWeight(G4double prob) {
  if (prob < 0.) prob = G4CMPConfigManager::GetGenCharges();

  // If prob=0., random throw always fails, never divides by zero
  return ((prob==1.) ? 1. : (G4UniformRand()<prob) ? 1./prob : 0.);
}

// Get current track from event and track managers

G4Track* G4CMP::GetCurrentTrack() {
  return G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack();
}

// Get touchable from current track

const G4VTouchable* G4CMP::GetCurrentTouchable() {
  G4Track* track = GetCurrentTrack();
  return track ? track->GetTouchable() : 0;
}

// Copy information from current step into data block]

void G4CMP::FillHit(const G4Step* step, G4CMPElectrodeHit* hit) {
  // Get information from the track
  G4Track* track     = step->GetTrack();
  G4int trackID      = track->GetTrackID();
  G4String name      = track->GetDefinition()->GetParticleName();
  G4double startE    = track->GetVertexKineticEnergy();
  G4double startTime = track->GetGlobalTime() - track->GetLocalTime();
  G4double finalTime = track->GetGlobalTime();
  G4double weight    = track->GetWeight();
  G4double edp       = step->GetNonIonizingEnergyDeposit();

  // Get start and end positions. Must use PreStepPoint to get correct
  // volume
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  G4ThreeVector startPosition = track->GetVertexPosition();
  G4ThreeVector finalPosition = postStepPoint->GetPosition();

  // Insert data into hit.
  hit->SetStartTime(startTime);
  hit->SetFinalTime(finalTime);
  hit->SetStartEnergy(startE);
  hit->SetEnergyDeposit(edp);
  hit->SetWeight(weight);
  hit->SetStartPosition(startPosition);
  hit->SetFinalPosition(finalPosition);
  hit->SetTrackID(trackID);
  hit->SetParticleName(name);
}


// Generate cos(theta) law for diffuse reflection, ensuring that computed
// vector is directed inward with respect to the surface normal.

G4ThreeVector
G4CMP::GetLambertianVector(const G4LatticePhysical* theLattice,
			   const G4ThreeVector& surfNorm, G4int mode) {
  const G4ThreeVector surfPoint = GetCurrentTrack()->GetPosition();
  return GetLambertianVector(theLattice, surfNorm, mode, surfPoint);
}

G4ThreeVector
G4CMP::GetLambertianVector(const G4LatticePhysical* theLattice,
			   const G4ThreeVector& surfNorm, G4int mode,
			   const G4ThreeVector& surfPoint) {
  G4ThreeVector reflectedKDir;
  const G4int maxTries = 1000;
  G4int nTries = 0;
  do {
    reflectedKDir = LambertReflection(surfNorm);
  } while (nTries++ < maxTries &&
           !PhononVelocityIsInward(theLattice, mode, reflectedKDir, surfNorm,
                                   surfPoint));

  return reflectedKDir;
}

G4ThreeVector G4CMP::LambertReflection(const G4ThreeVector& surfNorm) {
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1.0) / 2.0;

  G4ThreeVector refl = -surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}


// Check that phonon is properly directed from the volume surface
// waveVector and surfNorm need to be in global coordinates

G4bool G4CMP::PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int mode,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm) {
  const G4ThreeVector surfacePos = GetCurrentTrack()->GetPosition();
  return PhononVelocityIsInward(lattice, mode, waveVector, surfNorm, surfacePos);
}

G4bool G4CMP::PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int mode,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm,
                                     const G4ThreeVector& surfacePos) {
  // Get touchable for coordinate rotations
  const G4VTouchable* touchable = GetCurrentTouchable();

  if (!touchable) {
    G4Exception("G4CMP::PhononVelocityIsInward", "G4CMPUtils001",
		EventMustBeAborted, "Current track does not have valid touchable!");
    return false;
  }

  // MapKtoVDir requires local direction for the wavevector
  G4ThreeVector vDir = lattice->MapKtoVDir(mode, GetLocalDirection(touchable, waveVector));

  // Project a 1 nm step in the new direction, see if it
  // is still in the correct volume.
  G4ThreeVector localPos = GetLocalPosition(touchable, surfacePos);
  G4VSolid* solid = touchable->GetSolid();
  EInside trialStep = solid->Inside(localPos + 1*nm * vDir);

  // Compare group velocity and surface normal in global coordinates
  RotateToGlobalDirection(touchable, vDir);
  return (vDir.dot(surfNorm) < 0.0 && trialStep == kInside);
}


// Thermal distributions, useful for handling phonon thermalization

G4double G4CMP::MaxwellBoltzmannPDF(G4double temperature, G4double energy) {
  if (temperature <= 0.) return (energy == 0. ? 1. : 0.);

  const G4double kT = k_Boltzmann*temperature;

  // NOTE: coefficient usually has kT^-(3/2), but extra 1/kT makes units
  const G4double mbCoeff = 2. * sqrt(energy/(pi*kT));

  // This should be a true PDF, normalized to unit integral
  // NOTE: coefficient usually has kT^-(3/2), but extra 1/kT makes units
  return mbCoeff * exp(-energy/kT);
}

G4double G4CMP::ChooseThermalEnergy(G4double temperature) {
  // FIXME: With inverse CDF, we could do a simple direct throw
  // return G4CMP::MaxwellBoltzmannInvCDF(temperature, G4UniformRand());

  G4double kT = k_Boltzmann*temperature;
  G4double trialE;
  do {
    trialE = G4UniformRand() * 10.*kT;		// Uniform spread in E
  } while (!IsThermalized(temperature, trialE));

  return trialE;
}

G4double G4CMP::ChooseThermalEnergy(const G4LatticePhysical* lattice) {
  return lattice ? ChooseThermalEnergy(lattice->GetTemperature()) : 0.;
}


G4bool G4CMP::IsThermalized(G4double temperature, G4double energy) {
  return (G4UniformRand() < MaxwellBoltzmannPDF(temperature,energy));
}

G4bool G4CMP::IsThermalized(G4double energy) {
  return IsThermalized(G4CMPConfigManager::GetTemperature(), energy);
}

G4bool G4CMP::IsThermalized(const G4LatticePhysical* lattice, G4double energy) {
  return (lattice ? IsThermalized(lattice->GetTemperature(), energy)
	  : IsThermalized(energy) );		// Fall back to global temp.
}

G4bool G4CMP::IsThermalized(const G4Track* track) {
  if (!track) return false;
  return IsThermalized(G4CMP::GetLattice(*track), track->GetKineticEnergy());
}


// Search particle's processes for specified name

G4VProcess*
G4CMP::FindProcess(const G4ParticleDefinition* pd, const G4String& pname) {
  G4ProcessVector* pvec = pd->GetProcessManager()->GetProcessList();
  if (!pvec) return 0;

  for (size_t i=0; i<pvec->size(); i++) {
    G4VProcess* proc = (*pvec)[i];
    if (proc && proc->GetProcessName() == pname) return proc;
  }

  return 0;			// No match found
}


// Generate random index for shuffling secondaries

size_t G4CMP::RandomIndex(size_t n) {
  return (size_t)(n*G4UniformRand());
}


// Create debugging file with suffix or infix identifying worker thread

G4String G4CMP::DebuggingFileThread(const G4String& basefile) {
  if (!G4Threading::IsWorkerThread()) return basefile;	// No mods required

  G4String tid = "_G4WT" + std::to_string(G4Threading::G4GetThreadId());

  G4String tidfile = basefile;

  size_t lastdot = basefile.last('.');
  if (lastdot < basefile.length()) tidfile.insert(lastdot, tid);
  else tidfile += tid;

  return tidfile;
}
