/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPHitMerging.hh
/// \brief Definition of the G4CMPHitMerging factory class.  This class
///	is used by G4CMPSecondaryProduction, and may be used directly by
///     user applications, to consolidate adjacent "hits" (energy deposits)
///     along a track, for more efficient downsampling.
//
// $Id$
//
// 20220815  Michael Kelsey -- Extracted from G4CMPSecondaryProduction
// 20220821  G4CMP-308 -- Use new G4CMPStepInfo container instead of G4Step
// 20220826  For use with primary generator, need to pass in G4Event*
// 20220828  Add interface to process "left over" accumulators to primaries

#ifndef G4CMPHitMerging_hh
#define G4CMPHitMerging_hh 1

#include "G4CMPProcessUtils.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <vector>

class G4CMPEnergyPartition;
class G4Event;
class G4PrimaryParticle;
class G4Step;
class G4Track;
class G4VParticleChange;


class G4CMPHitMerging : public G4CMPProcessUtils {
public:
  G4CMPHitMerging();
  virtual ~G4CMPHitMerging();

  // Enable diagnostic messages
  void SetVerboseLevel(G4int vb=0) { verboseLevel = vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  // Configurable parameter to specify how to consolidate steps
  void SetCombiningStepLength(G4double val) { combiningStepLength = val; }
  G4bool GetCombiningStepLength() const { return combiningStepLength; }

  // Overload G4CMPProcessUtils function to fill energy parameters
  virtual void LoadDataForTrack(const G4Track* track);

  // Register "in process" G4Event*, for use when called from primary generator
  void ProcessEvent(const G4Event* primaryEvent=0);

  // Incorporate step into consolidated energy depost, generate secondaries
  // Return value indicate if new tracks are ready for use
  G4bool ProcessStep(const G4CMPStepInfo& stepData);
  G4bool ProcessStep(const G4Step& step);
  G4bool ProcessStep(const G4Step* step);

  // Transfer generated secondaries into process return object
  void FillOutput(G4VParticleChange* aParticleChange);

  // Transfer generated primaries into primary-track vertices
  void FillOutput(G4Event* primaryEvent, G4double time=0.);

  // Check for any unprocessed accumulators and convert them to primaries
  void FinishOutput(G4Event* primaryEvent);

protected:
  G4bool DoAddStep(const G4CMPStepInfo& step) const;	 // Accumulate step?
  G4bool HasEnergy(const G4CMPStepInfo& step) const;	 // Does step deposit?
  G4bool ReadyForOutput(const G4CMPStepInfo& step) const; // Ready to process?

  void PrepareOutput();		// Convert accumulator with partitioner

  // Process accumulator for track unconditionally to new primaries
  void FlushAccumulator(G4int trkID, G4Event* primaryEvent);

  // Create secondaries along the specified trajectory
  void GeneratePositions(size_t npos, const G4ThreeVector& start,
			 const G4ThreeVector& end);

  // Adjust position to be enforced inside current volume
  G4ThreeVector SurfaceClearance(const G4ThreeVector& pos);

private:
  G4int verboseLevel;			// Select diagnostic message output
  G4double combiningStepLength;		// Steps within which to accumulate
  G4bool readyForOutput;		// Flag if hit data ready for use

  // Collection of accumulators for individual tracks in event
  std::map<G4int, G4CMPStepAccumulator> trackAccum;
  G4CMPStepAccumulator* accumulator;	// Sums multiple steps along track
  G4int currentEventID;			// Remember event for clearing accums

  G4CMPEnergyPartition* partitioner;	// Creates secondary kinematics
  std::vector<G4ThreeVector> posSecs;	// Positions along trajectory
  std::vector<G4Track*> theSecs;		// Set of created secondaries
  std::vector<G4PrimaryParticle*> thePrims;	// Set of created primaries

  // No copying allowed
  G4CMPHitMerging(const G4CMPHitMerging& right);
  G4CMPHitMerging& operator=(const G4CMPHitMerging& right);
};

#endif	/* G4CMPHitMerging_hh */

