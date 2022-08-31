/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Converts phonon into QPs
#ifndef G4CMPQPDownconversion_h
#define G4CMPQPDownconversion_h 1

#include "G4VPhononProcess.hh"
#include "G4CMPBoundaryUtils.hh"
#include "G4CMPKaplanUtils.hh"

class G4CMPQPDownconversion : public G4VPhononProcess, public G4CMPKaplanUtils {
public:
  G4CMPQPDownconversion(const G4String& processName);
  ~G4CMPQPDownconversion();

  //virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
  //                                              G4double previousStepSize,
  //                                              G4ForceCondition* condition);
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

protected:
  //virtual G4double GetMeanFreePath(const G4Track& aTrack,
  //                                 G4double prevStepLength,
  //                                 G4ForceCondition* condition);

  G4bool IsSubgap(G4double energy) const { return energy < 2.*gapEnergy; }

private:
  // hide assignment operator as private
  G4CMPQPDownconversion(G4CMPQPDownconversion&);
  G4CMPQPDownconversion(G4CMPQPDownconversion&&);
  G4CMPQPDownconversion& operator=(const G4CMPQPDownconversion&);
  G4CMPQPDownconversion& operator=(const G4CMPQPDownconversion&&);

  void MakeQPSecondaries(const G4Track&);



private:
  G4int verboseLevel;		// For diagnostic messages

  G4MaterialPropertiesTable* filmProperties;
  G4double filmThickness;	// Quantities extracted from properties table
  G4double gapEnergy;		// Bandgap energy (delta)
  G4double subgapAbsorption;	// Probability to absorb energy below bandgap
  G4double vSound;		// Speed of sound in film

  mutable std::ofstream output;		// Diagnostic output under G4CMP_DEBUG
};

#endif	/* G4CMPPhononBoundaryProcess_h */
