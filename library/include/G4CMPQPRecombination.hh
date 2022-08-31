/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPQPRecombination_h
#define G4CMPQPRecombination_h 1

#include "globals.hh"
#include "G4CMPVProcess.hh"
#include "G4CMPKaplanUtils.hh"

class G4CMPQPRecombination : public G4CMPVProcess, public G4CMPKaplanUtils {
public:
  G4CMPQPRecombination();
  virtual ~G4CMPQPRecombination();

  //Using Rate models
  //using G4CMPVProcess::UseRateModel;		// Avoid function hiding
  //void UseRateModel(G4double model);		// Non-virtual to use in ctor

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  virtual bool IsApplicable(const G4ParticleDefinition&);

  G4CMPQPRecombination(G4CMPQPRecombination&) = delete;
  G4CMPQPRecombination(G4CMPQPRecombination&&) = delete;
  G4CMPQPRecombination& operator=(const G4CMPQPRecombination&) = delete;
  G4CMPQPRecombination& operator=(const G4CMPQPRecombination&&) = delete;

  static G4double GetMeanFreePath(const G4ParticleDefinition* pd);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:

};

#endif
