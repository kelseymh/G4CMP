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
#include "G4QP.hh"
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

G4bool G4CMP::IsQP(const G4Track& track) {
  return IsQP(track.GetParticleDefinition());
}

G4bool G4CMP::IsQP(const G4Track* track) {
  return (track!=0 && IsQP(*track));
}

G4bool G4CMP::IsQP(const G4ParticleDefinition& pd) {
  return IsQP(&pd);
}

G4bool G4CMP::IsQP(const G4ParticleDefinition* pd) {
  return (pd == G4QP::Definition());
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
			   const G4ThreeVector& surfNorm, G4int mode,
			   const G4ThreeVector& incidentVDir) {
  const G4ThreeVector surfPoint = GetCurrentTrack()->GetPosition();
  return GetLambertianVector(theLattice, surfNorm, mode, incidentVDir,
			     surfPoint);
}

// Note here that the surfNorm here can be a "naive" surface norm calculated
// from just considering the surface plane. The degeneracy between "inward"
// and "outward" is broken with the call to the generalizedsurfacenormal,
// which, when dotted into the incoming momentum, is always negative.

G4ThreeVector
G4CMP::GetLambertianVector(const G4LatticePhysical* theLattice,
			   const G4ThreeVector& surfNorm, G4int mode,
			   const G4ThreeVector& incidentVDir,
			   const G4ThreeVector& surfPoint) {

  G4ThreeVector reflectedKDir;
  const G4int maxTries = 10000;
  G4int nTries = 0;

  //Due to a need to generalize this lambertian vector finding, we need to
  //make sure that the incident momentum and the surface norm are in the same
  //directions. Otherwise, take surfNorm to be -1*surfNorm. If we don't
  //do this here, we can end up launching phonons in a lambertian way
  //but in the wrong direction relative to the surface, which may sometimes
  //be into vaccum.
  G4ThreeVector generalizedSurfaceNorm =
    G4CMP::GetGeneralizedSurfaceNormal(surfNorm,incidentVDir);

  
  do {
    reflectedKDir = LambertReflection(generalizedSurfaceNorm);
  } while (nTries++ < maxTries &&
           !PhononVelocityIsInward(theLattice, mode, reflectedKDir,
				   generalizedSurfaceNorm,
                                   surfPoint));

  //If we exceed our max tries, then set things to 0 so we know to kill the
  //track
  if (nTries >= maxTries) {
    //G4cout << "nTries >= maxTries in GetLambertianVector." << G4endl;
    reflectedKDir = G4ThreeVector(0,0,0);
  }
  
  return reflectedKDir;
}

//Now modified to recognize that the surface norm is "generalized," i.e. it is
//always pointing in the same direction as the incident velocity/momentum.
//This is opposite what existed before, where we could always assume we were
//reflecting "inward" with respect to an outward-facing surface normal.

G4ThreeVector G4CMP::LambertReflection(const G4ThreeVector& surfNorm) {
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1.0) / 2.0;

  G4ThreeVector refl = -1*surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}


// Check that phonon is properly directed from the volume surface
// waveVector and surfNorm need to be in global coordinates.
// A new assumption here is that the surfNorm that is passed in is "generalized,"
// i.e. that it points in the same direction as to the incoming v-dir (not
// the incoming k-vector). It therefore is negative when dotted into the outgoing
// v-dir (if reflection is occurring). So if you're using these, make sure that
// the surfNorm that you pass in is pointing *away from* the volume that your phonon
// is coming from. Otherwise you may confuse yourself.

G4bool G4CMP::PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int mode,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm) {
  const G4ThreeVector surfacePos = GetCurrentTrack()->GetPosition();
  return PhononVelocityIsInward(lattice, mode, waveVector, surfNorm,
				surfacePos);
}

// See above comment about the generalized surface norm being what is passed
// in here.

G4bool G4CMP::PhononVelocityIsInward(const G4LatticePhysical* lattice,
                                     G4int mode,
                                     const G4ThreeVector& waveVector,
                                     const G4ThreeVector& surfNorm,
                                     const G4ThreeVector& surfacePos) {
  // Get touchable for coordinate rotations
  const G4VTouchable* touchable = GetCurrentTouchable();

  if (!touchable) {
    G4Exception("G4CMP::PhononVelocityIsInward", "G4CMPUtils001",
		EventMustBeAborted,
		"Current track does not have valid touchable!");
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

  //G4cout << "In PVII: refl. wavevect: " << waveVector << ", vDir: " <<
  //  vDir << ", surfNorm: " << surfNorm << ", trialStep: " << trialStep <<
  //  ", kIn: " << kInside << G4endl;
  //G4cout << "In PVII: cur. touchable vol: " <<
  //  touchable->GetVolume()->GetName() << ", localPos: " << localPos <<
  //  ", trialStepPos: " << localPos + 1*nm * vDir << G4endl;
  
  //In all discernible cases, the (generalized) surface norm passed into this
  //function should now be pointing in the direction identical to the incident
  //velocity. If the new velocity dotted into this generalized surface norm
  //is negative (and the trial step is inside the incident volume), then we've
  //succeeded. Otherwise, return false. REL changed to < 0.0 8/25/25
  return (vDir.dot(surfNorm) < 0.0 && trialStep == kInside);
}


// Check that the phonon is properly directed outward from the surface it is
// impinging on. This is a bit distinct from the PhononVelocityIsInward
// function because in that one we get to assume that the transformations can
// all be done in the same touchable volume. In this, we cannot -- we must
// pass the next-volume's touchable to this function so that we can use it
// in the rotations from global to local and the mapping of K to vDir. You'll
// also need to make sure that the lattice is the *new* lattice (i.e. the one
// that the phonon is leaving into). As with PhononVelocityIsInward, this takes
// a "generalized" surfaceNorm, which always points in the same direction as
// the impinging phonon's velocity (not necessarily its k-vector). So make
// sure that's true as well. All of these input vectors should be in a global
// frame.

G4bool G4CMP::PhononVelocityIsOutward(const G4LatticePhysical* lattice,
				      G4int mode,
				      const G4ThreeVector& waveVector,
				      const G4ThreeVector& surfNorm,
				      const G4VTouchable * nextVolTouchable ){
  const G4ThreeVector surfacePos = GetCurrentTrack()->GetPosition();
  return PhononVelocityIsOutward(lattice, mode, waveVector, surfNorm,
				 nextVolTouchable, surfacePos);
}

// See above comments on details for surfNorm and lattice requirements.
G4bool G4CMP::PhononVelocityIsOutward(const G4LatticePhysical* lattice,
				      G4int mode,
				      const G4ThreeVector& waveVector,
				      const G4ThreeVector& surfNorm, 
				      const G4VTouchable * nextVolTouchable,
				      const G4ThreeVector& /*surfacePos*/) {
  // Get touchable for coordinate rotations
  const G4VTouchable* touchable = nextVolTouchable;
  
  if (!touchable) {
    G4Exception("G4CMP::PhononVelocityIsOutward", "G4CMPUtils000",
		EventMustBeAborted,
		"Current track does not have valid touchable!");
    return false;
  }

  // MapKtoVDir requires local direction for the wavevector
  G4ThreeVector vDir = lattice->MapKtoVDir(mode, GetLocalDirection(touchable, waveVector));
  
  // Compare group velocity and surface normal in global coordinates
  RotateToGlobalDirection(touchable, vDir);

  //G4cout << "In PhononVelocityIsOutward, vDir rotated to global direction is:
  //" << vDir << G4endl;
  
  //In all discernible cases, the (generalized) surface norm passed into this
  //function should now be pointing in the direction identical to the incident
  //velocity. If the new velocity dotted into this generalized surface norm is
  //positive, then we've succeeded and the phonon is "outbound". Otherwise,
  //return false. Note that we don't use a trialStep logic here because if this
  //is leaving a volume, it may either leave into a sibling volume (in which
  //case kInside == false, a mother volume (kInside == false), or a daughter
  //volume (kInside == true), so there's not a "hard and fast" logic to this
  //as there is in the "directed inward/reflection" case
  
  //G4cout << "In PhononVelocityIsOutward, Returning: " << (vDir.dot(surfNorm))
  //<< " with vDir: " << vDir << " and surfNorm: " << surfNorm << G4endl;
  return (vDir.dot(surfNorm) > 0.0); //REL changed 8/25/2025
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
