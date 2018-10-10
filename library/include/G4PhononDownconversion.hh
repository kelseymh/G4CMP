/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononDownconversion.hh
/// \brief Definition of the G4PhononDownconversion class
//
// $Id$
//
// 20170805  Replace GetMeanFreePath() with scattering-rate model

#ifndef G4PhononDownconversion_h
#define G4PhononDownconversion_h 1

#include "G4VPhononProcess.hh"

class G4CMPAnharmonicDecay;

class G4PhononDownconversion : public G4VPhononProcess {
public:
  G4PhononDownconversion(const G4String& processName ="phononDownconversion");
  virtual ~G4PhononDownconversion();
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

  // Only longitudinal phonons decay (L -> L' T, or L -> T T')
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond) {
    return G4CMPVProcess::GetMeanFreePath(trk, prevstep, cond);
  }


private:
  G4CMPAnharmonicDecay* anharmonicDecay;
  // hide assignment operator as private
  G4PhononDownconversion(G4PhononDownconversion&);
  G4PhononDownconversion& operator=(const G4PhononDownconversion& right);

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4PhononDownconversion_h */
