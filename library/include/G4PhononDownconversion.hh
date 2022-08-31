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
// 20181010  J. Singh -- Move functionality to G4CMPAnharmonicDecay.
// 20181011  M. Kelsey -- Add LoadDataForTrack() to initialize decay utility.
// 20201109  Drop G4CMP_DEBUG protection here, to avoid client rebuilding

#ifndef G4PhononDownconversion_h
#define G4PhononDownconversion_h 1

#include "G4VPhononProcess.hh"

class G4CMPAnharmonicDecay;


class G4PhononDownconversion : public G4VPhononProcess {
public:
  G4PhononDownconversion(const G4String& processName ="phononDownconversion");
  virtual ~G4PhononDownconversion();

  // Pass verbosity through to decay utility
  virtual void SetVerboseLevel(G4int vb);

  // Only longitudinal phonons decay (L -> L' T, or L -> T T')
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  // Configure for current track including AnharmonicDecay utility
  virtual void LoadDataForTrack(const G4Track* track);

  // Perform downconversion using AnharmonicDecay utility
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );

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
};

#endif	/* G4PhononDownconversion_h */
