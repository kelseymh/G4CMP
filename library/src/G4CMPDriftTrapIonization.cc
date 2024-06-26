/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20200331  G4CMP-196: Added impact ionization process
// 20200426  G4CMP-196: Change name to TrapIonization, specify beam and trap
//		particle types
// 20200604  G4CMP-208: Comment out unused function arguments

#include "G4CMPDriftTrapIonization.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <vector>


// For constructor, return abbreviation for particles, traps

namespace {
  const char* pdLetter(const G4ParticleDefinition* pd) {
    return (G4CMP::IsElectron(pd) ? "e" : G4CMP::IsHole(pd) ? "h" : "x");
  }

  const char* trapLetter(const G4ParticleDefinition* pd) {
    return (G4CMP::IsElectron(pd) ? "D" : G4CMP::IsHole(pd) ? "A" : "X");
  }

  G4CMPProcessSubType subtype(const G4ParticleDefinition* pd) {
    return (G4CMP::IsElectron(pd) ? fDTrapIonization :
	    G4CMP::IsHole(pd) ? fATrapIonization : fG4CMPProcess);
  }
}

// Constructor and destructor

G4CMPDriftTrapIonization::
G4CMPDriftTrapIonization(G4ParticleDefinition* impactPD,
			 G4ParticleDefinition* trapPD, const G4String &name)
  : G4CMPVDriftProcess(pdLetter(impactPD)+(trapLetter(trapPD)+name),
		       subtype(trapPD)), impactType(impactPD),
    trapType(trapPD) {;}

G4CMPDriftTrapIonization::~G4CMPDriftTrapIonization() {;}


// Compute MFP for specific beam and trap types

G4double G4CMPDriftTrapIonization::
GetMeanFreePath(const G4ParticleDefinition* impactPD,
		const G4ParticleDefinition* trapPD) {
  if (G4CMP::IsElectron(impactPD)) {
    return (G4CMP::IsElectron(trapPD) ? G4CMPConfigManager::GetEDTrapIonMFP() :
	    G4CMP::IsHole(trapPD) ? G4CMPConfigManager::GetEATrapIonMFP() :
	    DBL_MAX);
  } else if (G4CMP::IsHole(impactPD)) {
    return (G4CMP::IsElectron(trapPD) ? G4CMPConfigManager::GetHDTrapIonMFP() :
	    G4CMP::IsHole(trapPD) ? G4CMPConfigManager::GetHATrapIonMFP() :
	    DBL_MAX);
  } else {
    // FIXME:  Show we throw an exception here, or just return?
  }

  return DBL_MAX;	// Should never get here
}

G4double 
G4CMPDriftTrapIonization::GetMeanFreePath(const G4Track& aTrack, G4double,
					  G4ForceCondition* /*cond*/) {
  G4bool changedLattice = UpdateMeanFreePathForLatticeChangeover(aTrack);
  return GetMeanFreePath(impactType, trapType);
}


// Process action

G4VParticleChange* 
G4CMPDriftTrapIonization::PostStepDoIt(const G4Track& aTrack,
				       const G4Step& /*aStep*/) {
  aParticleChange.Initialize(aTrack);

  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepDoIt: "
           << aTrack.GetDefinition()->GetParticleName()
           << " impact ionization of a " << pdLetter(trapType)
	   << "-type impurity trap." << G4endl;
  }

  // NOTE: If secondary gets real energy/momentum transfer, we *MUST*
  //       do the proper kinematics and update aParticleChange.

  // Create secondary with minimal energy, assuming no momentum transfer
  G4Track* knockon = G4CMP::CreateSecondary(aTrack, trapType,
					    G4RandomDirection(), 1e-3*eV);
  aParticleChange.AddSecondary(knockon);

  ClearNumberOfInteractionLengthLeft();		// All processes should do this!
  return &aParticleChange;
}
