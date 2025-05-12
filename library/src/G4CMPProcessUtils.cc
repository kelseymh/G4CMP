/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

//
/// \file library/src/G4CMPProcessUtils.cc
/// \brief Implementation of the G4CMPProcessUtils class
///   Provides useful general functions to fetch and store lattice, access
///   and apply lattice parameters for phonons and charge carriers.
///
///   Use via multiple inheritance with concrete or base process classes
//
// $Id$
//
// 20140321  Move lattice-based placement transformations here, via Touchable
// 20140407  Add functions for phonon generation in Luke scattering
// 20140412  Add manual configuration options
// 20140509  Add ChoosePolarization() which uses DOS values from lattice
// 20141216  Set velocity "by hand" for secondary electrons
// 20150109  Use G4CMP_SET_ELECTRON_MASS to enable dynamic mass, velocity set
// 20150112  Add GetCurrentValley() function to get valley of current track,
//	     allow GetValley functions to treat holes, returning -1
// 20150309  Add Create*() functions which take position and energy arguments
//	     (for use with AlongStepDoIt() actions).
// 20150310  Fix CreateChargeCarrier to use momentum unit vector
// 20160610  Return regular (NOT Herring-Vogt) wave vector for electrons
// 20160809  BUG FIX:  th_phonon==0 is fine for computing energy.
// 20160825  Add assignment operators for cross-process configuration;
//	     move track identification functions to G4CMPUtils.
// 20160829  Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
// 20160906  Make GetSurfaceNormal() const.
// 20161004  Add new ChangeValley() function to avoid null selection
// 20170525  Drop explicit copy constructors; let compiler do the work
// 20170602  Local track identification functions apply to current track only
// 20170620  Drop local caching of transforms; call through to G4CMPUtils.
// 20170621  Drop local initialization of TrackInfo; StackingAction only
// 20170624  Improve initialization from track, use Navigator to infer volume
// 20201124  Change argument name in MakeGlobalRecoil() to 'krecoil' (track)
// 20201223  Add FindNearestValley() function to align electron momentum.
// 20210318  In LoadDataForTrack, kill a bad track, not the whole event.
// 20211001  "Reverse" Get*VelocityVector() and Get*Momentum() functions, so
//	     that G4Track::GetMomentum() is used even for electrons, and
//	     velocity is calculated from that.  Use internal vector buffer.
// 20211002  FindNearestValley() implementation moved to G4CMPGeometryUtils.
// 20211003  Add track touchable as data member, to create if needed
// 20230524  Expand GetCurrentTouchable() to create one for new tracks
// 20230807  Multiplied ChargeCarrierTimeStep by mach to get the correct
//		Luke scattering rate
// 20230808  Added derivation of ChargeCarrierTimeStep to explain bug fix.
// 20230831  Remove modifications to ChargeCarrierTimeStep(), they seem to
//		cause zero-length and NaN steps.
// 20240303  Add local currentTouchable pointer for non-tracking situations.
// 20240402  Drop FindTouchable() function.  Set currentTouchable internally
//		not available from track, and delete it at end of track.
// 20250124  Add FillParticleChange() to update phonon wavevector and Vg.
// 20250129  Rotate Vg in FillParticleChange() to global coordinates.
// 20250423  Add FillParticleChange() to update phonon position and touchable.
// 20250505  Update local time for phonon displacement in FillParticleChange.
// 20250508  Fix local and global coordinate system for phonon wavevectors.
// 20250512  Use tempvec2 for Vg in LoadDataForTrack to improve performance.

#include "G4CMPProcessUtils.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPParticleChangeForPhonon.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4RandomDirection.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"
#include "G4GeometryTolerance.hh"
#include "G4VSolid.hh"
#include <set>


// Constructor and destructor

G4CMPProcessUtils::G4CMPProcessUtils()
  : theLattice(nullptr), currentTrack(nullptr), currentVolume(nullptr),
    currentTouchable(nullptr) {;}

G4CMPProcessUtils::~G4CMPProcessUtils() {;}


// Initialization for current track

void G4CMPProcessUtils::LoadDataForTrack(const G4Track* track) {
  // WARNING!  This assumes track starts and ends in one single volume!
  SetCurrentTrack(track);
  SetLattice(track);

  if (!theLattice) {
    G4String msg = ("No lattice found for volume "
		    + track->GetVolume()->GetName()
		    + "; track will be killed");
    G4Exception("G4CMPProcessUtils::LoadDataForTrack", "Utils001",
		JustWarning, msg);

    // Force track to be killed immediately
    const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);

    return;	// No lattice, no special actions possible
  }

  // Sanity check -- track should already have kinematics container
  if (!G4CMP::HasTrackInfo(track)) {
    G4Exception("G4CMPProcessUtils::LoadDataForTrack", "Utils002",
		JustWarning, "No auxiliary info found for track");
    
    G4CMP::AttachTrackInfo(track);
  }

  // Transfer phonon wavevector into momentum direction for this step
  if (IsPhonon()) {
    G4CMPPhononTrackInfo* trackInfo =
      G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(*track);

    // Set momentum direction using already provided wavevector
    tempvec = trackInfo->k();
    RotateToLocalDirection(tempvec);

    const G4ParticleDefinition* pd = track->GetParticleDefinition();
    G4Track* tmp_track = const_cast<G4Track*>(track);

    // Set the momentum direction with the group velocity
    tempvec2 = theLattice->MapKtoVDir(G4PhononPolarization::Get(pd), tempvec);
    RotateToGlobalDirection(tempvec2);
    tmp_track->SetMomentumDirection(tempvec2);
  }
}


// Identify current track type to simplify some conditionals

G4bool G4CMPProcessUtils::IsPhonon() const {
  return G4CMP::IsPhonon(currentTrack);
}

G4bool G4CMPProcessUtils::IsElectron() const {
  return G4CMP::IsElectron(currentTrack);
}

G4bool G4CMPProcessUtils::IsHole() const {
  return G4CMP::IsHole(currentTrack);
}

G4bool G4CMPProcessUtils::IsChargeCarrier() const {
  return G4CMP::IsChargeCarrier(currentTrack);
}


// Cache volume associated with tracking, using Navigator if necessary

void G4CMPProcessUtils::SetCurrentTrack(const G4Track* track) {
  currentTrack = track;
  currentTouchable = nullptr;
  currentVolume = track ? track->GetVolume() : nullptr;

  if (!track) return;		// Avoid unnecessry work

  if (!currentVolume) {		// Primary tracks may not have volumes yet
    currentVolume = G4CMP::GetVolumeAtPoint(track->GetPosition());
    currentTouchable = G4CMP::CreateTouchableAtPoint(track->GetPosition());
    deleteTouchable = true;	// Avoid memory leak at end of track
  }
}

void G4CMPProcessUtils::SetLattice(const G4Track* track) {
  theLattice = track ? G4CMP::GetLattice(*track) : nullptr;
}


// Fetch lattice for specified volume, use in subsequent steps

void G4CMPProcessUtils::FindLattice(const G4VPhysicalVolume* volume) {
  currentVolume = volume;

  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(volume);

  if (!theLattice) {
    G4cerr << "WARNING: No lattice for volume " << volume->GetName() << G4endl;
  }
}


// Fill ParticleChange wavevector and group velocity for a given wavevector
// Wavevector is expected to be in the global coordinate frame
void G4CMPProcessUtils::FillParticleChange(G4ParticleChange& particleChange,
            const G4Track& track, const G4ThreeVector& wavevector) const {
  // Get phonon mode from track
  G4int mode = GetPolarization(track);

  // Get Vg from global wavevector
  G4ThreeVector vDir = theLattice->MapKtoVDir(mode, GetLocalDirection(wavevector));
  G4double v = theLattice->MapKtoV(mode, GetLocalDirection(wavevector));

  // Update trackInfo and particleChange
  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(track);
  trackInfo->SetWaveVector(wavevector);
  particleChange.ProposeVelocity(v);
  RotateToGlobalDirection(vDir);
  particleChange.ProposeMomentumDirection(vDir);
}

void G4CMPProcessUtils::FillParticleChange(G4CMPParticleChangeForPhonon& particleChange,
  const G4Step& step, const G4ThreeVector& position) const {
    // Update the local time to account for displacement
    G4double delta_t = (position - *particleChange.GetPosition()).mag() / particleChange.GetVelocity();
    G4StepPoint* postStep = step.GetPostStepPoint();
    particleChange.ProposeLocalTime(postStep->GetLocalTime() + delta_t);

    // Update position, touchable, and step status
    particleChange.ProposePosition(position);
    particleChange.ProposeTouchableHandle(step.GetPreStepPoint()->GetTouchableHandle());
    step.GetPostStepPoint()->SetStepStatus(fPostStepDoItProc);
}


// Delete current configuration before new track starts

void G4CMPProcessUtils::ReleaseTrack() {
  currentTrack = nullptr;
  currentVolume = nullptr;
  theLattice = nullptr;

  ClearTouchable();
}


// Register touchable owned by client code

void G4CMPProcessUtils::SetTouchable(const G4VTouchable* touch) {
  ClearTouchable();
  currentTouchable = touch;
  deleteTouchable = false;	// Client code retains ownership
}

void G4CMPProcessUtils::ClearTouchable() const {
  if (deleteTouchable) delete currentTouchable;
  currentTouchable = nullptr;
  deleteTouchable = false;
}


// Convert between local and global coordinates in currentTrack frame

G4ThreeVector
G4CMPProcessUtils::GetLocalDirection(const G4ThreeVector& dir) const {
  return G4CMP::GetLocalDirection(GetCurrentTouchable(),dir);
}

G4ThreeVector
G4CMPProcessUtils::GetLocalPosition(const G4ThreeVector& pos) const {
  return G4CMP::GetLocalPosition(GetCurrentTouchable(),pos);
}

void G4CMPProcessUtils::RotateToLocalDirection(G4ThreeVector& dir) const {
  return G4CMP::RotateToLocalDirection(GetCurrentTouchable(),dir);
}

void G4CMPProcessUtils::RotateToLocalPosition(G4ThreeVector& pos) const {
  return G4CMP::RotateToLocalPosition(GetCurrentTouchable(),pos);
}

// Convert between local and global coordinates in currentTrack frame

G4ThreeVector 
G4CMPProcessUtils::GetGlobalDirection(const G4ThreeVector& dir) const {
  return G4CMP::GetGlobalDirection(GetCurrentTouchable(),dir);
}

G4ThreeVector
G4CMPProcessUtils::GetGlobalPosition(const G4ThreeVector& pos) const {
  return G4CMP::GetGlobalPosition(GetCurrentTouchable(),pos);
}

void G4CMPProcessUtils::RotateToGlobalDirection(G4ThreeVector& dir) const {
  return G4CMP::RotateToGlobalDirection(GetCurrentTouchable(),dir);
}

void G4CMPProcessUtils::RotateToGlobalPosition(G4ThreeVector& pos) const {
  return G4CMP::RotateToGlobalPosition(GetCurrentTouchable(),pos);
}


// Access track position and momentum in local coordinates

G4ThreeVector G4CMPProcessUtils::GetLocalPosition(const G4Track& track) const {
  return GetLocalPosition(track.GetPosition());
}

void G4CMPProcessUtils::GetLocalPosition(const G4Track& track,
					 G4double pos[3]) const {
  tempvec = GetLocalPosition(track);
  pos[0] = tempvec.x();
  pos[1] = tempvec.y();
  pos[2] = tempvec.z();
}

G4ThreeVector G4CMPProcessUtils::GetLocalMomentum(const G4Track& track) const {
  return GetLocalDirection(track.GetMomentum());
}

void G4CMPProcessUtils::GetLocalMomentum(const G4Track& track, 
					 G4double mom[3]) const {
  tempvec = GetLocalMomentum(track);
  mom[0] = tempvec.x();
  mom[1] = tempvec.y();
  mom[2] = tempvec.z();
}

G4ThreeVector 
G4CMPProcessUtils::GetLocalVelocityVector(const G4Track& track) const {
  G4ThreeVector vel;

  if (G4CMP::IsElectron(track)) {
    vel = theLattice->MapPtoV_el(GetValleyIndex(track),GetLocalMomentum(track));
  } else {
    vel = track.CalculateVelocity() * track.GetMomentumDirection();
    RotateToLocalDirection(vel);
  }

  return vel;
}

void G4CMPProcessUtils::GetLocalVelocityVector(const G4Track &track,
                                               G4double vel[]) const {
  tempvec = GetLocalVelocityVector(track);
  vel[0] = tempvec.x();
  vel[1] = tempvec.y();
  vel[2] = tempvec.z();
}

G4ThreeVector G4CMPProcessUtils::GetLocalWaveVector(const G4Track& track) const {
  if (G4CMP::IsChargeCarrier(track)) {
    return GetLocalMomentum(track) / hbarc;
  } else if (G4CMP::IsPhonon(track)) {
    return GetLocalDirection(G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(track)->k());
  } else {
    G4Exception("G4CMPProcessUtils::GetLocalWaveVector", "DriftProcess002",
                EventMustBeAborted, "Unknown charge carrier");
    return G4ThreeVector();
  }
}

// Access track position and momentum in global coordinates
G4ThreeVector 
G4CMPProcessUtils::GetGlobalPosition(const G4Track& track) const {
  return track.GetPosition();
}

void G4CMPProcessUtils::GetGlobalPosition(const G4Track& track,
           G4double pos[3]) const {
  tempvec = GetGlobalPosition(track);
  pos[0] = tempvec.x();
  pos[1] = tempvec.y();
  pos[2] = tempvec.z();
}

G4ThreeVector 
G4CMPProcessUtils::GetGlobalMomentum(const G4Track& track) const {
  return track.GetMomentum();
}

void G4CMPProcessUtils::GetGlobalMomentum(const G4Track& track,
					  G4double mom[3]) const {
  tempvec = GetGlobalMomentum(track);
  mom[0] = tempvec.x();
  mom[1] = tempvec.y();
  mom[2] = tempvec.z();
}

G4ThreeVector 
G4CMPProcessUtils::GetGlobalVelocityVector(const G4Track& track) const {
  tempvec = GetLocalVelocityVector(track);
  RotateToGlobalDirection(tempvec);
  return tempvec;
}

void G4CMPProcessUtils::GetGlobalVelocityVector(const G4Track &track, G4double vel[]) const {
  tempvec = GetGlobalVelocityVector(track);
  vel[0] = tempvec.x();
  vel[1] = tempvec.y();
  vel[2] = tempvec.z();
}

G4double G4CMPProcessUtils::CalculateVelocity(const G4Track& track) const {
  if (IsElectron()) {
    G4ThreeVector ptrk = GetLocalMomentum(track);
    return GetLattice()->MapPtoV_el(GetValleyIndex(track), ptrk).mag();
  } else {
    return track.CalculateVelocity();
  }
}

G4double G4CMPProcessUtils::GetKineticEnergy(const G4Track &track) const {
  if (IsElectron()) {
    G4ThreeVector ptrk = GetLocalMomentum(track);
    return GetLattice()->MapPtoEkin(GetValleyIndex(track), ptrk);
  } else {  
    return track.GetKineticEnergy();
  }
}

// Return particle type for currently active track [set in LoadDataForTrack()]

const G4ParticleDefinition* G4CMPProcessUtils::GetCurrentParticle() const {
  return (currentTrack ? currentTrack->GetParticleDefinition() : 0);
}


// Return touchable for currently active track for transforms

const G4VTouchable* G4CMPProcessUtils::GetCurrentTouchable() const {
  if (!currentTrack) return currentTouchable;

  const G4VTouchable* touch = currentTrack->GetTouchable();
  if (touch) return touch;

  // Create a local touchable, to delete at end of track
  ClearTouchable();
  currentTouchable = G4CMP::CreateTouchableAtPoint(currentTrack->GetPosition());
  deleteTouchable = true;

  return currentTouchable;
}


// Access phonon particle-type/polarization indices

G4int G4CMPProcessUtils::GetPolarization(const G4Track& track) const {
  return G4PhononPolarization::Get(track.GetParticleDefinition());
}


// Generate random polarization from density of states

G4int G4CMPProcessUtils::ChoosePhononPolarization() const {
  return G4CMP::ChoosePhononPolarization(theLattice->GetLDOS(),
                                         theLattice->GetSTDOS(),
                                         theLattice->GetFTDOS());
}


// Convert K wave vector to track momentum

void G4CMPProcessUtils::MakeGlobalRecoil(G4ThreeVector& krecoil) const {
  // Convert recoil wave vector to momentum in local frame 
  if (IsElectron()) {
    krecoil = theLattice->MapKtoP(GetValleyIndex(GetCurrentTrack()),krecoil);
  } else if (IsHole()) {
    krecoil *= hbarc;
  } else {
    G4Exception("G4CMPProcessUtils::MakeGlobalPhonon", "DriftProcess006",
                EventMustBeAborted, "Unknown charge carrier");
  }

  RotateToGlobalDirection(krecoil);
}


// Generate direction angle for phonons in Luke scattering

G4double G4CMPProcessUtils::MakePhononTheta(G4double k, G4double ks) const {
  G4double u = G4UniformRand();
  G4double v = ks/k;
  G4double base = (1-u) * (1 - 3*v + 3*v*v - v*v*v); 	// (1-u)*(1-v)^3
  if (base < 0.0) return 0;
  
  G4double operand = v + pow(base, 1.0/3.0);   
  if (operand > 1.0) operand=1.0;
  
  return acos(operand);
}

// Compute energy of phonon in Luke Scattering

G4double G4CMPProcessUtils::MakePhononEnergy(G4double k, G4double ks,
					     G4double th_phonon) const {
  return MakePhononEnergy(2.*(k*cos(th_phonon)-ks));
}

G4double G4CMPProcessUtils::MakePhononEnergy(G4double q) const {
  return q * theLattice->GetSoundSpeed() * hbar_Planck;
}

// Compute direction angle for recoiling charge carrier

G4double G4CMPProcessUtils::MakeRecoilTheta(G4double k, G4double ks,
					    G4double th_phonon) const {
  if (th_phonon == 0.) return 0.;		// Avoid unnecessary work

  G4double kctks = k*cos(th_phonon) - ks;

  return acos( (k*k - 2*ks*kctks - 2*kctks*kctks)
	       / (k * sqrt(k*k - 4*ks*kctks)) );
}


// Generate random valley for charge carrier

G4int G4CMPProcessUtils::ChangeValley(G4int valley) const {
  // generate random valley offset (up to N-1)
  G4int nv = theLattice->NumberOfValleys();
  G4int dv = (G4int)(G4UniformRand()*(nv-1))+1;

  // Apply offset to change input to new value
  return (valley+dv) % nv;
}


// Access electron propagation direction/index

G4int G4CMPProcessUtils::GetValleyIndex(const G4Track& track) const {
  return G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(track)->ValleyIndex();
}

const G4RotationMatrix& 
G4CMPProcessUtils::GetValley(const G4Track& track) const {
  G4int iv = GetValleyIndex(track);
  return (iv>=0 ? theLattice->GetValley(iv) : G4RotationMatrix::IDENTITY);
}

// Find valley which aligns most closely with _local_ direction vector
G4int G4CMPProcessUtils::FindNearestValley(const G4Track& track) const {
  return (G4CMP::IsElectron(track) ? FindNearestValley(GetLocalMomentum(track))
	  : -1);
}

// NOTE:  Direction vector must be passed in _local_ coordinate system
G4int G4CMPProcessUtils::FindNearestValley(const G4ThreeVector& dir) const {
  return G4CMP::FindNearestValley(theLattice, dir);
}


// Compute characteristic time step for charge carrier
// Rate depends on Mach number (v/vsound = k/ksound) and scattering length l0
// 1/Tau = velLong/(3*l0) * kmag/ksound * (1 - ksound/kmag)^3
// Tau   = (3*l0)/velLong * 1/mach * (1 - 1/mach)^-3
//       = (3*l0)/velLong * 1/mach * [(mach-1)/mach]^-3
//       = (3*l0)/velLong * mach^3/mach * (mach-1)^-3
//       = (3*l0)/velLong * mach^2 / (mach-1)^3

G4double 
G4CMPProcessUtils::ChargeCarrierTimeStep(G4double mach, G4double l0) const {
  const G4double velLong = theLattice->GetSoundSpeed();

  const G4double tstep = 3.*l0/velLong;
  return (mach<1.) ? tstep : tstep*mach/((mach-1)*(mach-1)*(mach-1));
  // NOTE: Above numerator should be tstep*mach*mach, but causes problems
}
