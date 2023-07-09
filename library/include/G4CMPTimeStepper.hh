/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
//	     Add wrapper function to compute individual time steps
// 20140418  Remove local valley transforms; use lattice functions
// 20150112  Drop redundant IsApplicable (identical to base version)
// 20170905  Cache Luke and IV rate models in local LoadDataFromTrack()
// 20170908  Drop "time step" functions, use rate models as estimators
// 20190906  Provide functions to externally set rate models
// 20200331  C. Stanford (G4CMP-195): Added charge trapping
// 20200331  G4CMP-196: Added impact ionization mean free path
// 20200426  G4CMP-196: Change "impact" name to "trapIon"
// 20220730  G4CMP-301: Drop trapping processes, as they have built-in MFPs,
//		don't need TimeStepper for energy-dependent calculation.
// 20230527  Drop competitive MFP calculation; use this process to enforce
//	       a maximum allowed step length, to support mass recalculation.

#ifndef G4CMPTimeStepper_h
#define G4CMPTimeStepper_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"

class G4CMPTimeStepper : public G4CMPVDriftProcess {
public:
  G4CMPTimeStepper();
  virtual ~G4CMPTimeStepper();

  // No random throw here: MFP and GPIL are fixed lengths
  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double prevStepSize,
				       G4ForceCondition* condition) {
    return GetMeanFreePath(aTrack, prevStepSize, condition);
  }

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
					  const G4Step& aStep);

protected:  
  virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*);

private:
  //hide assignment operator
  G4CMPTimeStepper(G4CMPTimeStepper&);
  G4CMPTimeStepper& operator=(const G4CMPTimeStepper& right);
};

#endif
