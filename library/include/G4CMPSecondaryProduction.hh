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

#ifndef G4CMPSecondaryProduction_hh
#define G4CMPSecondaryProduction_hh 1

#include "G4VContinuousProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <utility>

class G4Track;
class G4Step;
class G4VParticleChange;
class G4ParticleDefinition;


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
  virtual G4double GetContinuousStepLimit(const G4Track& aTrack,
					  G4double  previousStepSize,
					  G4double  currentMinimumStep,
					  G4double& currentSafety);

protected:
  void AddPhonons(const G4Step& stepData);
  
  void AddChargeCarriers(const G4Step& stepData);
  
  G4int GenerateEnergyPositions(const G4Step& stepData, G4double Etotal,
				G4double yield, G4double sigma);

private:
  // NOTE:  These are temporary local params -- will move to lattice config
  G4double ionizationEnergy;		// Get this from G4Material
  G4double yieldPerPhonon;
  G4double yieldPerChargePair;
  G4double sigmaPerPhonon;
  G4double sigmaPerChargePair;

  typedef std::pair<G4double, G4ThreeVector> EPosPair;
  std::vector<EPosPair> energyPosList;		// Buffer for secondaries
  
  // No copying allowed
  G4CMPSecondaryProduction(const G4CMPSecondaryProduction& right);
  G4CMPSecondaryProduction& operator=(const G4CMPSecondaryProduction& right);
};

#endif	/* G4CMPSecondaryProduction_hh */
