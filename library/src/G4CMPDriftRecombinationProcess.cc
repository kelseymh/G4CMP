/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPDriftRecombinationProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"

G4CMPDriftRecombinationProcess::G4CMPDriftRecombinationProcess(
                                                       const G4String &name,
                                                       G4CMPProcessSubType type)
                          : G4VRestProcess(name, fPhonon), G4CMPProcessUtils() {
  verboseLevel = G4CMPConfigManager::GetVerboseLevel();
  SetProcessSubType(type);
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4bool G4CMPDriftRecombinationProcess::IsApplicable(
                                              const G4ParticleDefinition& aPD) {
  return (&aPD==G4CMPDriftElectron::Definition() ||
          &aPD==G4CMPDriftHole::Definition());
}

void G4CMPDriftRecombinationProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  LoadDataForTrack(track);
}

void G4CMPDriftRecombinationProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  ReleaseTrack();
}

G4double G4CMPDriftRecombinationProcess::AtRestGetPhysicalInteractionLength(
                                                  const G4Track& aTrack,
                                                  G4ForceCondition* condition) {
  return GetMeanLifeTime(aTrack, condition);
}

G4double G4CMPDriftRecombinationProcess::GetMeanLifeTime(
                                              const G4Track& /*aTrack*/,
                                              G4ForceCondition* /*condition*/) {
  // Charge carriers should immediately recombine when KE = 0.
  return 0.;
}

G4VParticleChange* G4CMPDriftRecombinationProcess::AtRestDoIt(
                                                      const G4Track& aTrack,
                                                      const G4Step& aStep) {
  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftRecombinationProcess::AtRestDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " reabsorbed by the lattice."
           << G4endl;
  }

  aParticleChange.Initialize(aTrack);

  // FIXME: Each charge carrier is independent, so it only gives back 0.5 times
  // the band gap. Really electrons and holes should recombine, killing both
  // tracks and giving back the band gap. Maybe there is a better way?

  // FIXME: What does the recombo phonon distribution look like?

  G4double weight = G4CMP::ChoosePhononWeight();
  if (weight > 0.) {
    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,
                                   G4RandomDirection(),
                                   0.5 * theLattice->GetBandGapEnergy());
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);
  }
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  return &aParticleChange;
}
