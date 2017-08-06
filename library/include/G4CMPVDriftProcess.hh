/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVDriftProcess.hh
/// \brief Definition of the G4CMPVDriftProcess base class
//
// $Id$
//
// 20140312  Introduce multiple inheritance from G4CMPProcessUtils
// 20140325  Move time-step calculation here from TimeStepper and LukeScat
// 20140331  Add required subtype code to constructor
// 20140902  Add new kinematics function which takes energy as input
// 20141231  Add function to enforce minimum step length (fraction of L0)
// 20150112  Rename SetNewKinematics to FillParticleChange for clarity
// 20150601  Inherit from new G4CMPVProcess
// 20170806  Move ChargeCarrierTimeStep() to G4CMPProcessUtils

#ifndef G4CMPVDriftProcess_h
#define G4CMPVDriftProcess_h 1

#include "G4CMPVProcess.hh"
#include "G4CMPProcessSubType.hh"
#include "G4ThreeVector.hh"


class G4CMPVDriftProcess : public G4CMPVProcess {
public:
  G4CMPVDriftProcess(const G4String& processName, G4CMPProcessSubType stype);
  virtual ~G4CMPVDriftProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

  // Overload base version to set a minimum step size, avoiding "stuck" tracks
  virtual G4double
  PostStepGetPhysicalInteractionLength(const G4Track& track,
				       G4double previousStepSize,
				       G4ForceCondition* condition);

protected:
  // Fill ParticleChange energy and mass for charge carrier of given momentum
  void FillParticleChange(G4int ivalley, const G4ThreeVector& p);

  // Fill ParticleChange energy and mass for charge carrier of given energy
  void FillParticleChange(G4int ivalley, G4double Ekin,
              const G4ThreeVector& v);

private:
  // hide assignment operators as private 
  G4CMPVDriftProcess(G4CMPVDriftProcess&);
  G4CMPVDriftProcess& operator=(const G4CMPVDriftProcess& right);
};

#endif	/* G4CMPVDriftProcess_h */
