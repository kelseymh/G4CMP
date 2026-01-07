/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPairBreakingRate.hh
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
///
/// This is the rate class for the QP diffusion time stepper process.
/// Please see G4CMPQPDiffusionTimeStepperProcess.hh for more info.
/// The .cc file corresponding to this is (for now) where the rate
/// for this process is hardcoded.
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPQPDiffusionTimeStepperRate_hh
#define G4CMPQPDiffusionTimeStepperRate_hh 1

#include "G4CMPVScatteringRate.hh"

class G4CMPQPDiffusionTimeStepperRate : public G4CMPVScatteringRate {
public:
  G4CMPQPDiffusionTimeStepperRate() : G4CMPVScatteringRate("qpDiffusionTimeStepper") {;}
  virtual ~G4CMPQPDiffusionTimeStepperRate() {;}
  
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical* theLat);

private:
};

#endif	/* G4CMPQPDiffusionTimeStepperRate_hh */
