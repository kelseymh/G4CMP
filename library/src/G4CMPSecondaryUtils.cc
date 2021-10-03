/// \file library/include/G4CMPSecondaryUtils.hh
/// \brief Free standing helper functions for creating secondary tracks
///
//
// $Id$
//
// 20161115 Initial commit - R. Agnese
// 20170620 M. Kelsey -- Replace PV arg with Touchable, for transforms
// 20170629 M. Kelsey -- Add volume name to "no lattice" error messages.
// 20170721 M. Kelsey -- Check volume in AdjustSecondaryPosition.
// 20170815 M. Kelsey -- Move AdjustSecondaryPosition to GeometryUtils
// 20170928 M. Kelsey -- Replace "polarization" with "mode"
// 20210518 M. Kelsey -- Protect new secondaries from production cuts
// 20211001 M. Kelsey -- Collapse layered CreateChargeCarrier functions

#include "G4CMPSecondaryUtils.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"


// Generic function to create both phonon and charge carrier secondaries

G4Track* G4CMP::CreateSecondary(const G4Track& track,
                                G4ParticleDefinition* pd,
                                const G4ThreeVector& waveVec,
                                G4double energy) {
  const G4VTouchable* touch = track.GetTouchable();

  if (G4CMP::IsPhonon(pd)) {
    return CreatePhonon(touch, G4PhononPolarization::Get(pd),
                        waveVec, energy, track.GetGlobalTime(), track.GetPosition());
  }

  if (G4CMP::IsChargeCarrier(pd)) {
    // Assign valley in which momentum is best aligned
    G4int valley = (G4CMP::IsElectron(pd)
		    ? G4CMP::FindNearestValley(touch, waveVec) : -1);

    return CreateChargeCarrier(touch, G4int(pd->GetPDGCharge()/eplus),
                               valley, energy, track.GetGlobalTime(),
                               waveVec.unit(), track.GetPosition());
  }

  G4Exception("G4CMP::CreateSecondary", "Secondary001", EventMustBeAborted,
              ("Particle Definition "+pd->GetParticleName()
	       +" does not match a G4CMP particle.").c_str());
  return nullptr;
}


// Create phonon with specified mode and wave vector

G4Track* G4CMP::CreatePhonon(const G4VTouchable* touch, G4int mode,
                             const G4ThreeVector& waveVec,
                             G4double energy, G4double time,
			     const G4ThreeVector& pos) {
  G4LatticePhysical* lat = G4CMP::GetLattice(touch);
  if (!lat) {
    G4Exception("G4CMP::CreatePhonon", "Secondary002", EventMustBeAborted,
                ("No lattice for volume "+touch->GetVolume()->GetName()).c_str());
    return nullptr;
  }

  if (mode == G4PhononPolarization::UNKNOWN) {		// Choose random value
    mode = ChoosePhononPolarization(lat);
  }

  G4ThreeVector vgroup = lat->MapKtoVDir(mode, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cerr << "WARNING: vgroup not a unit vector: " << vgroup
	   << " length " << vgroup.mag() << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(mode);

  // Secondaries are (usually) created at the current track coordinates
  RotateToGlobalDirection(touch, vgroup);

  auto sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
                         time, G4CMP::ApplySurfaceClearance(touch, pos));
  sec->SetGoodForTrackingFlag(true);	// Protect against production cuts

  // Store wavevector in auxiliary info for track
  AttachTrackInfo(sec, GetGlobalDirection(touch, waveVec));

  sec->SetVelocity(lat->MapKtoV(mode, waveVec));
  sec->UseGivenVelocity(true);

  return sec;
}

G4Track* G4CMP::CreateChargeCarrier(const G4VTouchable* touch, G4int charge,
                                    G4int valley, G4double Ekin, G4double time,
                                    const G4ThreeVector& pdir,
                                    const G4ThreeVector& pos) {
  if (charge != 1 && charge != -1) {		// Sanity check
    G4Exception("G4CMP::CreateChargeCarrier","Secondary004",EventMustBeAborted,
                "Invalid charge for charge carrier.");
    return nullptr;
  }

  G4LatticePhysical* lat = G4CMP::GetLattice(touch);
  if (!lat) {
    G4Exception("G4CMP::CreateChargeCarrier","Secondary003",EventMustBeAborted,
                ("No lattice for volume "+touch->GetVolume()->GetName()).c_str());
    return nullptr;
  }

  // If electron wasn't given a valley, use best alignment with momentum
  if (charge < 0 && valley == -1) valley = G4CMP::FindNearestValley(touch,pdir);

  G4ParticleDefinition* theCarrier = nullptr;
  G4double carrierMass=0.;

  if (charge > 0) {
    theCarrier  = G4CMPDriftHole::Definition();
    carrierMass = lat->GetHoleMass();
  } else {
    theCarrier  = G4CMPDriftElectron::Definition();
    G4ThreeVector plocal = G4CMP::GetLocalDirection(touch,pdir);
    carrierMass = lat->GetElectronEffectiveMass(valley, plocal); 
  }

  // NOTE:  G4CMP uses true mass unts: convert MeV/c^2 to MeV for Geant4
  G4DynamicParticle* secDP = new G4DynamicParticle(theCarrier, pdir.unit(),
						   Ekin, carrierMass*c_squared);

  G4Track* sec = new G4Track(secDP, time, G4CMP::ApplySurfaceClearance(touch, pos));
  sec->SetGoodForTrackingFlag(true);	// Protect against production cuts

  // Store wavevector in auxiliary info for track
  G4CMP::AttachTrackInfo(sec, valley);

  // Temporary warning about hole valleys
  if (valley >= 0 && G4CMP::IsHole(theCarrier)) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary010", JustWarning,
		"Hole has been assigned a valley index.");
  }

  return sec;
}
