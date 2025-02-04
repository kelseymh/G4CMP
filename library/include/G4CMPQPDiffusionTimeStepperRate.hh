/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPairBreakingRate.hh
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
//

#ifndef G4CMPQPDiffusionTimeStepperRate_hh
#define G4CMPQPDiffusionTimeStepperRate_hh 1

#include "G4CMPVScatteringRate.hh"
#include "G4CMPSCUtils.hh"
#include <vector>
#include <map>

class G4CMPQPDiffusionTimeStepperRate : public G4CMPVScatteringRate {
					   
public:
  G4CMPQPDiffusionTimeStepperRate() : G4CMPVScatteringRate("QPDiffusionTimeStepper") {;}
  virtual ~G4CMPQPDiffusionTimeStepperRate() {;}
  
  virtual G4double Rate(const G4Track& aTrack) const;
  virtual void UpdateLookupTable(const G4LatticePhysical * theLat);
  

private:
  
  
};

#endif	/* G4CMPQPDiffusionTimeStepperRate_hh */
