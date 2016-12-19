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
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

G4Track* G4CMP::CreateSecondary(const G4Track& track,
                                G4ParticleDefinition* pd,
                                const G4ThreeVector& waveVec,
                                G4double energy) {
  if (G4CMP::IsPhonon(pd)) {
    return CreatePhonon(track.GetVolume(), G4PhononPolarization::Get(pd),
                        waveVec, energy, track.GetGlobalTime(), track.GetPosition());
  }

  if (G4CMP::IsChargeCarrier(pd)) {
    const G4LatticePhysical* lat = GetTrackInfo<G4CMPDriftTrackInfo>(track)->Lattice();
    return CreateChargeCarrier(track.GetVolume(), G4int(pd->GetPDGCharge()/eplus),
                               ChooseValley(lat), energy, track.GetGlobalTime(),
                               waveVec, track.GetPosition());
  }

  G4Exception("G4CMP::CreateTrack", "Secondary001", EventMustBeAborted,
              "Particle Definition does not match a G4CMP particle.");
  return nullptr;
}

G4Track* G4CMP::CreatePhonon(const G4VPhysicalVolume* vol, G4int polarization,
                             const G4ThreeVector& waveVec,
                             G4double energy, G4double time, const G4ThreeVector& pos) {
  G4ThreadLocalStatic auto latMan = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = latMan->GetLattice(vol);
  if (!lat) {
    G4Exception("G4CMP::CreatePhonon", "Secondary002", EventMustBeAborted,
                "No lattice for volume.");
  }

  if (polarization == G4PhononPolarization::UNKNOWN) {		// Choose value
    polarization = ChoosePhononPolarization(lat);
  }

  G4ThreeVector vgroup = lat->MapKtoVDir(polarization, waveVec);
  if (std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cerr << "WARNING: vgroup not a unit vector: " << vgroup
     << " length " << vgroup.mag() << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(polarization);

  // Secondaries are (usually) created at the current track coordinates
  RotateToGlobalDirection(vol, vgroup);
  G4ThreeVector secPos = AdjustSecondaryPosition(vol, pos);

  auto sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
                         time, secPos);

  // Store wavevector in auxiliary info for track
  auto secInfo = new G4CMPPhononTrackInfo(lat, GetGlobalDirection(
                                            vol, waveVec));
  AttachTrackInfo(*sec, secInfo);

  sec->SetVelocity(lat->MapKtoV(polarization, waveVec));
  sec->UseGivenVelocity(true);

  return sec;
}

G4Track* G4CMP::CreateChargeCarrier(const G4VPhysicalVolume* vol, G4int charge,
                                    G4int valley, G4double Ekin, G4double time,
                                    const G4ThreeVector& pdir,
                                    const G4ThreeVector& pos) {
  G4ThreadLocalStatic auto latMan = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = latMan->GetLattice(vol);
  if (!lat) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary003", EventMustBeAborted,
                "No lattice for volume.");
  }

  if (charge != 1 && charge != -1) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary004", EventMustBeAborted,
                "Invalid charge for charge carrier.");
  }

  if (charge == 1) { // Hole
    G4double pmag = std::sqrt(2. * Ekin * lat->GetHoleMass());
    return CreateChargeCarrier(vol, charge, valley, time, pmag * pdir, pos);
  } else if (charge == -1) { // Electron
    G4double k_HVmag = std::sqrt(2. * Ekin * lat->GetElectronMass()) / hbar_Planck;
    G4ThreeVector k_HVdir = lat->MapV_elToK_HV(valley, pdir).unit();
    G4ThreeVector p = lat->MapK_HVtoP(valley, k_HVmag * k_HVdir);
    return CreateChargeCarrier(vol, charge, valley, time, p, pos);
  }

  G4Exception("G4CMP::CreateChargeCarrier", "Secondary005", EventMustBeAborted,
              "Couldn't recognize charge carrier.");
  return nullptr;
}

G4Track* G4CMP::CreateChargeCarrier(const G4VPhysicalVolume* vol, G4int charge,
                                    G4int valley, G4double time,
                                    const G4ThreeVector& p,
                                    const G4ThreeVector& pos) {
  G4ThreadLocalStatic auto latMan = G4LatticeManager::GetLatticeManager();
  G4LatticePhysical* lat = latMan->GetLattice(vol);
  if (!lat) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary006", EventMustBeAborted,
                "No lattice for volume.");
  }

  if (charge != 1 && charge != -1) {
    G4Exception("G4CMP::CreateChargeCarrier", "Secondary007", EventMustBeAborted,
                "Invalid charge for charge carrier.");
  }

  G4ParticleDefinition* theCarrier = nullptr;
  G4double carrierMass=0., carrierEnergy=0.;

  G4ThreeVector v_unit(0.);
  if (charge == 1) {
    theCarrier    = G4CMPDriftHole::Definition();
    carrierMass   = lat->GetHoleMass();
    carrierEnergy = 0.5 * p.mag2() / carrierMass;	// Non-relativistic
    v_unit = p.unit();
  } else {
    theCarrier    = G4CMPDriftElectron::Definition();
    carrierMass   = lat->GetElectronMass();
    G4ThreeVector p_local = G4CMP::GetLocalDirection(vol, p);
    G4ThreeVector v_local = lat->MapPtoV_el(valley, p_local);
    G4CMP::RotateToGlobalDirection(vol, v_local); // v_local is now actually global
    carrierEnergy = 0.5 * carrierMass * v_local.mag2();// Non-relativistic
    v_unit = v_local.unit();
  }

  auto secDP = new G4DynamicParticle(theCarrier, v_unit, carrierEnergy, carrierMass);

  G4ThreeVector secPos = G4CMP::AdjustSecondaryPosition(vol, pos);
  auto sec = new G4Track(secDP, time, secPos);

  // Store wavevector in auxiliary info for track
  auto trackInfo = new G4CMPDriftTrackInfo(lat, valley);
  G4CMP::AttachTrackInfo(*sec, trackInfo);

  return sec;
}

G4ThreeVector G4CMP::AdjustSecondaryPosition(const G4VPhysicalVolume* vol,
                                             G4ThreeVector pos) {
  // Take a copy because we would've had to make a copy at some point anyway.
  // If the step is near a boundary, create the secondary in the initial volume
  G4Navigator* nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // Safety is the distance to a boundary
  G4double safety = nav->ComputeSafety(pos);
  // Tolerance is the error in deciding which volume a track is in
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  // If the distance to an edge is within error, we might accidentally get placed
  // in the next volume. Instead, let's scoot a bit away from the edge.
  if (safety <= kCarTolerance) {
    G4ThreeVector norm = vol->GetLogicalVolume()->GetSolid()
                         ->SurfaceNormal(GetLocalPosition(vol, pos));
    RotateToGlobalDirection(vol, norm);
    pos += (safety - kCarTolerance * (1.001)) * norm;
  }

  return pos;
}