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

#ifndef G4CMPSecondaryProduction_hh
#define G4CMPSecondaryProduction_hh 1

#include "G4VContinuousProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4CMPEnergyPartition;
class G4DynamicParticle;
class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VParticleChange;


class G4CMPSecondaryProduction : public G4VContinuousProcess,
				 public G4CMPProcessUtils {
public:
  G4CMPSecondaryProduction();
  virtual ~G4CMPSecondaryProduction();

  // Applies to all charged, non-resonance particles
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  // Generate secondaries based on already-computed energy loss
  virtual G4VParticleChange* AlongStepDoIt(const G4Track& track,
					   const G4Step& stepData);

  // Overload G4CMPProcessUtils function to fill energy parameters
  virtual void LoadDataForTrack(const G4Track* track);

protected:
  // Calculate step limit for Along Step (not needed here)
  virtual G4double GetContinuousStepLimit(const G4Track&, G4double, G4double,
					  G4double&);

protected:
  void AddSecondaries(const G4Step& stepData);
  void GeneratePositions(const G4Step& stepData, size_t npos);

public:
  static size_t RandomIndex(size_t imax);	// Used to randomize secondaries

private:
  G4CMPEnergyPartition* partitioner;		// Creates secondary kinematics
  std::vector<G4Track*> theSecs;		// List of created secondaries
  std::vector<G4ThreeVector> posSecs;		// Positions along trajectory

  // No copying allowed
  G4CMPSecondaryProduction(const G4CMPSecondaryProduction& right);
  G4CMPSecondaryProduction& operator=(const G4CMPSecondaryProduction& right);
};

#endif	/* G4CMPSecondaryProduction_hh */
