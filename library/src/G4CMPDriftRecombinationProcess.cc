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
  : G4CMPVDriftProcess(name, type) {
}

G4VParticleChange* G4CMPDriftRecombinationProcess::PostStepDoIt(
                                                      const G4Track& aTrack,
                                                      const G4Step&) {
  aParticleChange.Initialize(aTrack);

  // If the particle has not come to rest, do nothing
  if (aTrack.GetTrackStatus() != fStopButAlive) {
    return &aParticleChange;
  }

  if (verboseLevel > 1) {
    G4cout << "G4CMPDriftRecombinationProcess::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " reabsorbed by the lattice."
           << G4endl;
  }

  // FIXME: Each charge carrier is independent, so it only gives back 0.5 times
  // the band gap. Really electrons and holes should recombine, killing both
  // tracks and giving back the band gap. Maybe there is a better way?

  // FIXME: What does the recombo phonon distribution look like?
  // For now we'll just divvy up the gap energy into n Debye energy phonons.
  G4double ePot = 0.5 * theLattice->GetBandGapEnergy();
  G4double eDeb = theLattice->GetDebyeEnergy();
  size_t n = std::ceil(ePot / eDeb);
  aParticleChange.SetNumberOfSecondaries(n);
  while (ePot > 0.) {
    G4double E = ePot > eDeb ? eDeb : ePot;
    ePot -= eDeb;
    G4Track* phonon = CreatePhonon(G4PhononPolarization::UNKNOWN,
                                   G4RandomDirection(), E);
    aParticleChange.AddSecondary(phonon);
  }

  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}

G4double G4CMPDriftRecombinationProcess::GetMeanFreePath(const G4Track&,
                                                         G4double,
                                                         G4ForceCondition* cond) {
  *cond = Forced;
  return DBL_MAX;
}

