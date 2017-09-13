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

#ifndef G4CMPTimeStepper_h
#define G4CMPTimeStepper_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"

class G4CMPVScatteringRate;


class G4CMPTimeStepper : public G4CMPVDriftProcess {
public:
  G4CMPTimeStepper();
  virtual ~G4CMPTimeStepper();

  // Initialize local pointers to Luke and IV scattering rate models
  virtual void LoadDataForTrack(const G4Track* aTrack);

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

  const G4double minStep;		// Minimum global step length 

  // Maximum rate for other processes, given track kinematics
  G4double MaxRate(const G4Track& aTrack) const;

  // Find distance in E-field to get to Luke threshold (Vsound)
  G4double StepToLuke(const G4Track& aTrack) const;

  // Energy increase due to field at current location
  G4double EnergyStep(const G4Track& aTrack, G4double step) const;

  // Compute track kinematics applying field acceleration
  void AdjustKinematics(const G4Track& aTrack, G4double deltaE);
  void CopyTrack(const G4Track& aTrack);
  G4Track* tempTrack;

  // Get scattering rates for other charge-carrier processes
  void ReportRates(const G4Track& aTrack);
  const G4CMPVScatteringRate* lukeRate;		// Rate models for current track
  const G4CMPVScatteringRate* ivRate;

private:
  //hide assignment operator
  G4CMPTimeStepper(G4CMPTimeStepper&);
  G4CMPTimeStepper& operator=(const G4CMPTimeStepper& right);
};

#endif
