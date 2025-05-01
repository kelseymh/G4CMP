/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPIVSwitchModel.hh
/// \brief Simple but non-physical implementation of intervalley scattering
///        of electrons during charge transport.  Reassigns valley index
///        registered for track, moving the momentum direction to preserve
///	   the angle with respect to the valley axis (conserves energy).
//
// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140418  Drop valley transforms; get from lattice
// 20170802  Remove MFP calculation; use scattering-rate model
// 20190704  Add selection of rate model by name, and material specific
// 20190906  For rate model selection, pass string by value
// 20190906  Push selected rate model back to G4CMPTimeStepper for consistency
// 20250430  Move PostStepDoIt() implementation from G4CMPInterValleyScattering

#ifndef G4CMPIVSwitchModel_h
#define G4CMPIVSwitchModel_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPIVSwitchModel : public G4CMPVPhysicsModel { 
public:
  G4CMPIVSwitchModel(G4CMPVProcess* process)
    : G4CMPVPhysicsModel("IVSwitch", process) {;}
  virtual ~G4CMPIVSwitchModel();

  // Do scattering action here
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

private:
  //hide assignment operator as private
  G4CMPIVSwitchModel(G4CMPIVSwitchModel&);
  G4CMPIVSwitchModel& operator=(const G4CMPIVSwitchModel& right);
};

#endif
