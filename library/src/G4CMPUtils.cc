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

#include "G4CMPUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4CMPBogoliubov.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"


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

G4bool G4CMP::IsQuasiparticle(const G4Track& track) {
  return IsQuasiparticle(track.GetParticleDefinition());
}

G4bool G4CMP::IsQuasiparticle(const G4Track* track) {
  return (track!=0 && IsQuasiparticle(*track));
}

G4bool G4CMP::IsQuasiparticle(const G4ParticleDefinition& pd) {
  return IsQuasiparticle(&pd);
}

G4bool G4CMP::IsQuasiparticle(const G4ParticleDefinition* pd) {
  return (pd == G4CMPBogoliubov::Definition());
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


// Generate cos(theta) law for diffuse reflection

G4ThreeVector G4CMP::LambertReflection(const G4ThreeVector& surfNorm) {
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1.0) / 2.0;

  G4ThreeVector refl = -surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}


// Check that phonon is properly directed from the volume surface

G4bool G4CMP::PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int mode,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm) {
  G4ThreeVector vDir = lattice->MapKtoVDir(mode, waveVector);
  return vDir.dot(surfNorm) < 0.0;
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
