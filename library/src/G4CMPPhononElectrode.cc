/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// class G4CMPPhononElectrode, subclass of G4CMPVElectrodePattern
//
// Class description:
//
// This provides an interface to G4CMPKaplanQP, to simulate the response
// of a phonon energy-absorbing electrode made of a thin-film superconductor.
//
// The associated border surface must have a G4MaterialPropertiesTable
// registered, containing at least all the following parameters (see also
// README.md):
//
// | Property Key        | Definition                   | Example value (Al) |
// |---------------------|------------------------------|--------------------|
// | filmThickness       | Thickness of film            | 600.*nm            |
// | gapEnergy           | Bandgap of film material     | 173.715e-6*eV      |
// | lowQPLimit          | Minimum bandgap multiple     | 3.                 |
// | phononLifetime      | Phonon lifetime at 2*bandgap | 242.*ps            |
// | phononLifetimeSlope | Lifetime vs. energy          | 0.29               |
// | vSound              | Speed of sound in film       | 3.26*km/s          |
//
// In addition, for sensors containing a "direct absorber" (such as a TES),
// The property key "subgapAbsorption" may be set with the probability to
// directly absorb phonons below 2*bandgap.
// 
// 20221006  M. Kelsey -- Adapted from SuperCDMS simulation version

#include "G4CMPPhononElectrode.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPKaplanQP.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"


// Constructor and destructor

G4CMPPhononElectrode::G4CMPPhononElectrode()
  : G4CMPVElectrodePattern(), kaplanQP(0) {;}

G4CMPPhononElectrode::~G4CMPPhononElectrode() {
  delete kaplanQP; kaplanQP=0;
}

// Assumes that user has configured a border surface only at sensor pads

G4bool G4CMPPhononElectrode::IsNearElectrode(const G4Step& /*step*/) const {
  return G4UniformRand() < GetMaterialProperty("filmAbsorption");
}


// Phonon gets killed upon breaking a Cooper pair. New phonons may
// be created and emitted back into the crystal

void G4CMPPhononElectrode::
AbsorbAtElectrode(const G4Track& track, const G4Step& step,
		  G4ParticleChange& particleChange) const {
  if (verboseLevel>1) {
    G4cout << "G4CMPPhononElectrode::AbsorbAtElectrode: Track "
	   << track.GetTrackID() << " absorbed at step "
	   << step.GetTrack()->GetCurrentStepNumber()
           << G4endl;
  }

  // Create KaplanQP simulator if not already available
  if (!kaplanQP) {
    G4MaterialPropertiesTable* ncTable =
      const_cast<G4MaterialPropertiesTable*>(theSurfaceTable);
    kaplanQP = new G4CMPKaplanQP(ncTable, verboseLevel);
  }

  // Transfer phonon energy into superconducting film
  G4double Ekin = GetKineticEnergy(track);

  // Track deposits some energy and may spawn new phonons (filled in KaplanQP)
  phononEnergies.clear();
  G4double EDep = kaplanQP->AbsorbPhonon(Ekin, phononEnergies);

  if (verboseLevel>1) {
    G4cout << " Incident Ekin " << Ekin << " absorbed " << EDep
	   << " emitted " << phononEnergies.size() << " secondaries" << G4endl;
  }

  if (EDep > 0. || phononEnergies.size() > 1) {
    ProcessAbsorption(track, step, EDep, particleChange);
  } else if (phononEnergies.size() == 1) {
    ProcessReflection(track, step, particleChange);
  } else {
    // Incident phonon is killed, no energy deposition
    particleChange.ProposeTrackStatus(fStopAndKill);
    particleChange.ProposeEnergy(0.);
  }
}


// Record energy deposition and re-emitted energies as secondary phonons

void G4CMPPhononElectrode::
ProcessAbsorption(const G4Track& track, const G4Step& step, G4double EDep,
		  G4ParticleChange& particleChange) const {
  // Incident phonon is killed, no energy deposition
  particleChange.ProposeTrackStatus(fStopAndKill);
  particleChange.ProposeEnergy(0.);

  particleChange.ProposeNonIonizingEnergyDeposit(EDep);

  // Secondaries are emitted with cos(theta) distribution inward
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(step);

  // Create secondaries for all of the generated phonon energies
  particleChange.SetNumberOfSecondaries(phononEnergies.size());

  G4double Ekin = GetKineticEnergy(track);
  G4ThreeVector k = GetLocalWaveVector(track);

  G4ThreeVector reflectedKDir;
  for (G4double E : phononEnergies) {
    G4double kmag = k.mag()*E/Ekin;	// Scale k vector by energy
    G4int pol = ChoosePhononPolarization();
    do {
      reflectedKDir = G4CMP::LambertReflection(surfNorm);
    } while (!G4CMP::PhononVelocityIsInward(theLattice, pol,
                                            kmag*reflectedKDir, surfNorm));

    G4Track* phonon = G4CMP::CreatePhonon(GetCurrentTouchable(),
					  pol, kmag*reflectedKDir,
					  E, track.GetGlobalTime(),
					  track.GetPosition());
    particleChange.AddSecondary(phonon);
  }	// for (E : ...)

  // Sanity check: secondaries' energy should equal assigned E
  if (verboseLevel>1) {
    G4double Esum = 0.;
    for (G4int i=0; i<particleChange.GetNumberOfSecondaries(); i++) {
      Esum += particleChange.GetSecondary(i)->GetKineticEnergy();
      G4cout << " secondary " << i << " GetKineticEnergy() "
	     << particleChange.GetSecondary(i)->GetKineticEnergy()/eV << " eV"
	     << " weight " << particleChange.GetSecondary(i)->GetWeight()
	     << G4endl;
    }

    if (fabs(Ekin-EDep-Esum) > 1e-6) {
      G4cerr << "ERROR: Energy non-conservation: Ekin-EDep " << Ekin-EDep
	     << " vs. sum of secondaries " << Esum << G4endl;
    }
  }
}


// Reflect phonon by changing direction randomly, no energy deposition

void G4CMPPhononElectrode::
ProcessReflection(const G4Track& track, const G4Step& step,
		  G4ParticleChange& particleChange) const {
  // Reflection is relative to inward surface normal
  G4ThreeVector surfNorm = G4CMP::GetSurfaceNormal(step);

  G4int pol = GetPolarization(track);

  G4ThreeVector reflectedKDir;
  do {
    reflectedKDir = G4CMP::LambertReflection(surfNorm);
  } while (!G4CMP::PhononVelocityIsInward(theLattice, pol,
					  reflectedKDir, surfNorm));

  if (verboseLevel>1)
    G4cout << " Phonon reflected from QET toward " << reflectedKDir << G4endl;

  particleChange.ProposeMomentumDirection(reflectedKDir);

  auto trackInfo = G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(track);
  trackInfo->IncrementReflectionCount();
}
