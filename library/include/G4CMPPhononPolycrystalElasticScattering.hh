/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhononPolycrystalElasticScattering.hh
/// \brief Definition of the G4CMPPhononPolycrystalElasticScattering class
///
/// This class defines the process for phonon elastic scattering on
/// effective grain boundaries in a polycrystalline substrate. A scale
/// of grain boundary is provided as a parameter, and scattering is
/// isotropic in k-space.
///
// $Id$
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPPhononPolycrystalElasticScattering_h
#define G4CMPPhononPolycrystalElasticScattering_h 1

#include "G4VPhononProcess.hh"

class G4CMPPhononPolycrystalElasticScattering : public G4VPhononProcess {
public:
  G4CMPPhononPolycrystalElasticScattering(const G4String& processName=
				       "phononPolycrystalElasticScattering");
  virtual ~G4CMPPhononPolycrystalElasticScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // Keep function here as call-back to avoid getting old toolkit version
  virtual G4double GetMeanFreePath(const G4Track& trk, G4double prevstep,
				   G4ForceCondition* cond);

private:
  // hide assignment operator as private 
  G4CMPPhononPolycrystalElasticScattering(G4CMPPhononPolycrystalElasticScattering&);
  G4CMPPhononPolycrystalElasticScattering& operator=
  (const G4CMPPhononPolycrystalElasticScattering& right);
};

#endif	/* G4CMPPhononPolycrystalElasticScattering_h */
