/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSecondaryProduction.hh
/// \brief Definition of the G4CMPSecondaryProduction process class.  This
///	class will be used to generate phonons and charge carriers as
///     secondaries, based on energy loss along step.
//
// $Id$
//
// 20150310  Michael Kelsey
// 20160825  Replace implementation with use of G4CMPEnergyPartition
// 20201207  Add flag to suspend parent track for secondary processing.
// 20210203  G4CMP-241 : Process must run after PostStepDoIt, not AlongStep.
// 20210303  G4CMP-243 : Consolidate nearby steps into one effective hit.
// 20210318  G4CMP-245 : Enforce clearance from crystal surfaces.
// 20210608  G4CMP-260 : Add function to identify steps with energy deposit.
// 20210610  G4CMP-262 : Handle step accumulation including track suspension,
//	       by keeping a map of accumulators by track ID
// 20220216  G4CMP-290 : Add start/end arguments to GeneratePositions().
// 20220815  G4CMP-308 : Factor step-accumulation procedures to HitMerging.

#ifndef G4CMPSecondaryProduction_hh
#define G4CMPSecondaryProduction_hh 1

#include "G4CMPVProcess.hh"

class G4CMPHitMerging;
class G4ParticleDefinition;
class G4Step;
class G4Track;


class G4CMPSecondaryProduction : public G4CMPVProcess {
public:
  G4CMPSecondaryProduction();
  virtual ~G4CMPSecondaryProduction();

  // Applies to all charged, non-resonance particles
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  // Generate secondaries based on already-computed energy loss
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
					   const G4Step& stepData);

  // Overload G4CMPProcessUtils function to fill energy parameters
  virtual void LoadDataForTrack(const G4Track* track);

  // Configurable flag to suspend parent track and process secondaries
  void ProcessSecondariesFirst(G4bool val) { secondariesFirst = val; }
  G4bool ProcessSecondariesFirst() const { return secondariesFirst; }

protected:
  G4double GetMeanFreePath(const G4Track&, G4double,
			   G4ForceCondition* condition);

private:
  G4CMPHitMerging* mergeHits;		// Utility to accumulate steps
  G4bool secondariesFirst;		// Process secondaries immediately

  // No copying allowed
  G4CMPSecondaryProduction(const G4CMPSecondaryProduction& right);
  G4CMPSecondaryProduction& operator=(const G4CMPSecondaryProduction& right);
};

#endif	/* G4CMPSecondaryProduction_hh */
