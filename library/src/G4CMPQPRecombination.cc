/***********************************************************************\
 *  * This software is licensed under the terms of the GNU General Public *
 *   * License version 3 or later. See G4CMP/LICENSE for the full license. *
 *   \***********************************************************************/


#include "G4CMPQPRecombination.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4RandomDirection.hh"
#include "G4Track.hh"
#include <vector>

G4CMPQPRecombination::G4CMPQPRecombination()
    : G4CMPVProcess("G4CMPQPRecombination",fQPRecombination)  {
    //UseRateModel(G4CMPConfigManager::GetQPRecombinationMFP());
}



G4CMPQPRecombination::~G4CMPQPRecombination() {;}


G4double 
G4CMPQPRecombination::GetMeanFreePath(const G4ParticleDefinition* pd) {
  return (G4CMP::IsQuasiparticle(pd) ? G4CMPConfigManager::GetQPRecombinationMFP()
	  : DBL_MAX);
}

G4double G4CMPQPRecombination::GetMeanFreePath(const G4Track&, G4double,
						    G4ForceCondition*) {
  return GetMeanFreePath(GetCurrentParticle());
}

/*
G4double 
G4CMPQPRecombination::GetMeanFreePath(const G4Track&, G4double,
						G4ForceCondition* cond) {
    *cond = Forced;
    return DBL_MAX;
}*/



G4bool 
G4CMPQPRecombination::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsQuasiparticle(&aPD);
}

G4VParticleChange* 
G4CMPQPRecombination::PostStepDoIt(const G4Track& aTrack,
					     const G4Step& aStep) {
    aParticleChange.Initialize(aTrack);

    //once QP has travelled, recombines with another QP

    G4double qpE = aTrack.GetKineticEnergy();
    G4double qpE2 = QPEnergyRand(2*theLattice->GetBandGapEnergy()); //Needs to be converted to a DOS calculation


    //Emits Phonon and Excess Energy is deposited into Lattice
    //PhononE = qp1E + qp2E - cooper-pairE
    G4double ePot = qpE+qpE2 - (2*theLattice->GetBandGapEnergy());
    //ePot*=0.5 ? As we only have 1 track?
    

    if (aParticleChange.GetNumberOfSecondaries() == 0) {	// Record energy release
      aParticleChange.ProposeNonIonizingEnergyDeposit(ePot);
    }

    aParticleChange.ProposeTrackStatus(fStopAndKill); //Kill Track

    ClearNumberOfInteractionLengthLeft();		// All processes should do this!

    return &aParticleChange;
}


