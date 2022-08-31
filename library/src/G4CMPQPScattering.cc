/***********************************************************************\
 *  * This software is licensed under the terms of the GNU General Public *
 *   * License version 3 or later. See G4CMP/LICENSE for the full license. *
 *   \***********************************************************************/


#include "G4CMPQPScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4Track.hh"
#include <vector>

G4CMPQPScattering::
G4CMPQPScattering()
    : G4CMPVProcess("G4CMPQPRecombination",fG4CMPProcess)  {
    UseRateModel(G4CMPConfigManager::GetQPRecombinationMFP());
}




G4CMPQPScattering::~G4CMPQPScattering() {;}

G4double 
G4CMPQPScattering::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
    *cond = Forced;
    return DBL_MAX;
}

G4bool 
G4CMPQPScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsQuasiparticle(&aPD);
}

G4VParticleChange* 
G4CMPQPScattering::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  //once QP has travelled, can scatter with emmision/absorb of phonon

  G4double qpE = aTrack.GetKineticEnergy();

  if(G4UniformRand() < 0.5) {//Should be diff number?

      //QP emits phonon
    G4double phononE = PhononEnergyRand(qpE);
    G4double qpENew = qpE-phononE; 

    if (aParticleChange.GetNumberOfSecondaries() == 0) {	// Record energy release
      aParticleChange.ProposeNonIonizingEnergyDeposit(phononE);
    }

    if (IsSubgap(qpENew)) { 
      aParticleChange.ProposeNonIonizingEnergyDeposit(qpENew);
      aParticleChange.ProposeTrackStatus(fStopAndKill); //Kill Track
    } else {
      aParticleChange.ProposeEnergy(qpENew);
    }
  } else {
      //absorbs Phonon
      //Choose Phonon from DOS
      G4double phononE = 0; //TBD
      
      //Change mom to random (no function for this?)
      aParticleChange.ProposeEnergy(qpE+phononE);      

    }

    ClearNumberOfInteractionLengthLeft();		// All processes should do this!
    return &aParticleChange;
}
