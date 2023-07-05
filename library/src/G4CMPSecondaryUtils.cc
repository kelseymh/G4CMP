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
// 20220907 G4CMP-316 -- Pass track into CreateXYZ() functions; do valley
//		selection for electrons in CreateChargeCarrier().

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


G4Track* G4CMP::CreateSecondary(const G4Track& track, G4ParticleDefinition* pd,
                                const G4ThreeVector& waveVec, G4double energy) {
  if (G4CMP::IsPhonon(pd)) {
    return CreatePhonon(track, G4PhononPolarization::Get(pd), waveVec,
                        energy, track.GetGlobalTime(), track.GetPosition());
  }

  if (G4CMP::IsChargeCarrier(pd)) {
    return CreateChargeCarrier(track, G4int(pd->GetPDGCharge()/eplus),
                               -1, energy, track.GetGlobalTime(),
                               waveVec, track.GetPosition());
  }

  G4Exception("G4CMP::CreateSecondary", "Secondary001", EventMustBeAborted,
              ("Particle Definition "+pd->GetParticleName()
	       +" does not match a G4CMP particle.").c_str());
  return nullptr;
}

G4Track* G4CMP::CreatePhonon(const G4Track& track, G4int mode,
			     const G4ThreeVector& waveVec, G4double energy,
			     G4double time, const G4ThreeVector& pos) {
  G4LatticePhysical* lat = G4CMP::GetLattice(track);
  if (!lat) {
    G4Exception("G4CMP::CreatePhonon", "Secondary002", EventMustBeAborted,
                ("No lattice for volume "+track.GetVolume()->GetName()).c_str());
    return nullptr;
  }

  if (mode == G4PhononPolarization::UNKNOWN) {		// Choose value
    mode = ChoosePhononPolarization(lat);
  }

  G4ThreeVector vgroup = lat->MapKtoVDir(mode, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cerr << "WARNING: vgroup not a unit vector: " << vgroup
     << " length " << vgroup.mag() << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(mode);

  // Secondaries are (usually) created at the current track coordinates
  const G4VTouchable* touch = track.GetTouchable();
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

// DEPRECATED: Version used by application code with output of G4CMPKaplanQP

G4Track* G4CMP::CreatePhonon(const G4VTouchable* touch, G4int mode,
			     const G4ThreeVector& waveVec, G4double energy,
			     G4double time, const G4ThreeVector& pos) {
  G4VPhysicalVolume* vol = touch->GetVolume();
  G4ThreadLocalStatic auto latMan = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = latMan->GetLattice(vol);
  if (!lat) {
    G4Exception("G4CMP::CreatePhonon", "Secondary002", EventMustBeAborted,
                ("No lattice for volume "+vol->GetName()).c_str());
    return nullptr;
  }

  if (mode == G4PhononPolarization::UNKNOWN) {		// Choose value
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


G4Track* G4CMP::CreateChargeCarrier(const G4Track& track, G4int charge,
                                    G4int valley, G4double Ekin, G4double time,
                                    const G4ThreeVector& pdir,
                                    const G4ThreeVector& pos) {

  G4LatticePhysical* lat = G4CMP::GetLattice(track);
  if (!lat) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary003", EventMustBeAborted,
                ("No lattice for volume "+track.GetVolume()->GetName()).c_str());
    return nullptr;
  }

  G4ThreeVector p = pdir;
  if (charge == 1) { 				// Hole
    p *= std::sqrt(2.*Ekin*lat->GetHoleMass());
    valley = -1;
  } else if (charge == -1) {			// Electron
    G4double k_HVmag = std::sqrt(2.*Ekin*lat->GetElectronMass()) / hbar_Planck;
    G4ThreeVector k_HVdir = lat->MapV_elToK_HV(valley, pdir).unit();
    p = lat->MapK_HVtoP(valley, k_HVmag * k_HVdir);
    if (valley<0) valley = ChooseValley(lat);
  } else {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary004", EventMustBeAborted,
                "Invalid charge for charge carrier.");
    return nullptr;
  }


  return CreateChargeCarrier(track, charge, valley, time, p, pos);
}

G4Track* G4CMP::CreateChargeCarrier(const G4Track& track, G4int charge,
                                    G4int valley, G4double time,
                                    const G4ThreeVector& p,
                                    const G4ThreeVector& pos) {
  G4LatticePhysical* lat = G4CMP::GetLattice(track);
  if (!lat) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary006", EventMustBeAborted,
                ("No lattice for volume "+track.GetVolume()->GetName()).c_str());
    return nullptr;
  }

  if (charge != 1 && charge != -1) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary007", EventMustBeAborted,
                "Invalid charge for charge carrier.");
    return nullptr;
  }

  const G4VTouchable* touch = track.GetTouchable();
  G4ParticleDefinition* theCarrier = nullptr;
  G4double carrierMass=0., carrierEnergy=0.;

  G4ThreeVector v_unit(0.);
  if (charge == 1) {
    theCarrier    = G4CMPDriftHole::Definition();
    carrierMass   = lat->GetHoleMass();
    carrierEnergy = 0.5 * p.mag2() / carrierMass;  // Non-relativistic
    v_unit = p.unit();
  } else {
    theCarrier    = G4CMPDriftElectron::Definition();
    carrierMass   = lat->GetElectronMass();
    G4ThreeVector p_local = G4CMP::GetLocalDirection(touch, p);
    G4ThreeVector v_local = lat->MapPtoV_el(valley, p_local);
    RotateToGlobalDirection(touch, v_local); // v_local is now actually global
    carrierEnergy = 0.5 * carrierMass * v_local.mag2();// Non-relativistic
    v_unit = v_local.unit();
    if (valley<0) valley = ChooseValley(lat);
  }

  // NOTE:  We use true mass unts: convert e.g. MeV/c^2 to MeV here
  auto secDP = new G4DynamicParticle(theCarrier, v_unit, carrierEnergy,
				     carrierMass*c_squared);

  auto sec = new G4Track(secDP, time, G4CMP::ApplySurfaceClearance(touch, pos));
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
