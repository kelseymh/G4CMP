/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPUtils.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"

G4CMPElectrodeSensitivity::G4CMPElectrodeSensitivity(G4String name)
  :G4VSensitiveDetector(name), hitsCollection(nullptr) {
  collectionName.insert("G4CMPElectrodeHit");
}

G4CMPElectrodeSensitivity::G4CMPElectrodeSensitivity(G4CMPElectrodeSensitivity&& in) :
  G4VSensitiveDetector(std::move(in)),
  hitsCollection(in.hitsCollection) {
}

G4CMPElectrodeSensitivity& G4CMPElectrodeSensitivity::operator=(G4CMPElectrodeSensitivity&& in) {
  // Base class members
  SensitiveDetectorName = std::move(in.SensitiveDetectorName);
  thePathName = std::move(in.thePathName);
  fullPathName = std::move(in.fullPathName);
  verboseLevel = in.verboseLevel;
  active = in.active;
  std::swap(ROgeometry, in.ROgeometry);
  std::swap(filter, in.filter);

  // Our members
  hitsCollection = in.hitsCollection;

  return *this;
}

void G4CMPElectrodeSensitivity::Initialize(G4HCofThisEvent* HCE) {
  hitsCollection = new G4CMPElectrodeHitsCollection(SensitiveDetectorName,
                                                    collectionName[0]);
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  HCE->AddHitsCollection(HCID, hitsCollection);
}

G4bool G4CMPElectrodeSensitivity::ProcessHits(G4Step* aStep,
                                              G4TouchableHistory* ROhist) {
  if (IsHit(aStep, ROhist)) {
    auto hit = new G4CMPElectrodeHit;
    G4CMP::FillHit(aStep, hit); // Mutates hit
    hitsCollection->insert(hit);
  }

  return true;
}

G4bool G4CMPElectrodeSensitivity::IsHit(const G4Step* step,
                                        const G4TouchableHistory*) const {
  /* Charge carriers do not deposit energy when they land on an electrode.
   * Phonons tracks are sometimes killed at the boundary in order to spawn new
   * phonon tracks. These tracks that are killed deposit no energy and should
   * not be picked up as hits.
   */
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition* particle = track->GetDefinition();

  G4bool isCharge = G4CMP::IsChargeCarrier(particle);
  G4bool isPhonon = G4CMP::IsPhonon(particle);

  G4bool deadAtBoundary = step->GetTrack()->GetTrackStatus() == fStopAndKill &&
                          postStepPoint->GetStepStatus() == fGeomBoundary;

  G4bool deposited = step->GetNonIonizingEnergyDeposit() > 0.;

  return deadAtBoundary && (isCharge || (isPhonon && deposited));
}
