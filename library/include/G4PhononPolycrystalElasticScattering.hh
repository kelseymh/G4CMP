/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononPolycrystalElasticScattering.hh
/// \brief Definition of the G4PhononPolycrystalElasticScattering class
//
// $Id$
//
// 20170805  Move GetMeanFreePath() to scattering-rate model

#ifndef G4PhononPolycrystalElasticScattering_h
#define G4PhononPolycrystalElasticScattering_h 1

#include "G4VPhononProcess.hh"

class G4PhononPolycrystalElasticScattering : public G4VPhononProcess {
public:
  G4PhononPolycrystalElasticScattering(const G4String& processName=
				       "phononPolycrystalElasticScattering");
  virtual ~G4PhononPolycrystalElasticScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond);

private:
  // hide assignment operator as private 
  G4PhononPolycrystalElasticScattering(G4PhononPolycrystalElasticScattering&);
  G4PhononPolycrystalElasticScattering& operator=
  (const G4PhononPolycrystalElasticScattering& right);
};

#endif	/* G4PhononPolycrystalElasticScattering_h */
