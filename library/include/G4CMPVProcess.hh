/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVProcess.hh
/// \brief Top-level base class for all G4CMP physics processes.  Only
///        discrete (post-step) processes are supported in G4CMP.  This
///	   base class provides access to all of the "ProcessUtils"
///	   functionaltiy via multiple inheritance.  Concrete processes
///	   don't need to do anything special.
//
// $Id$
//
// 20170601  New abstract base class for all G4CMP processes
// 20170802  Add registration of external scattering rate (MFP) model
// 20170905  Add accessors to get currentlty active scattering rate
// 20190906  Add function to initialize rate model after LoadDataForTrack
// 20250430  Add ability to register different physics models, which will
//	       be called by (new) base implementation of PostStepDoIt().

#ifndef G4CMPVProcess_h
#define G4CMPVProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPProcessUtils.hh"

class G4CMPVPhysicsModel;
class G4CMPVScatteringRate;


class G4CMPVProcess : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4CMPVProcess(const G4String& processName, G4CMPProcessSubType stype);
  virtual ~G4CMPVProcess();

  // Register utility class for computing scattering rate for MFP
  // NOTE:  Takes ownership of model for deletion
  // Subclasses MAY overload this to register a physics model
  virtual void UseRateModel(G4CMPVScatteringRate* model);
  const G4CMPVScatteringRate* GetRateModel() const { return rateModel; }
        G4CMPVScatteringRate* GetRateModel()       { return rateModel; }

  // Register specific physics model implementing kinematics
  void UsePhysicsModel(G4CMPVPhysicsModel* model);
  const G4CMPVPhysicsModel* GetPhysicsModel() const { return physicsModel; }
        G4CMPVPhysicsModel* GetPhysicsModel()       { return physicsModel; }

  // Initialize track/volume information (lattice, wavevector, etc.)
  // NOTE:  Subclasses must call back to these base implementations!
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

  // Calls through to physics model unconditionally
  // Subclasses MAY override, but if they do, they must not call back here
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  void ConfigureRateModel();		// Subclasses can call this directly

  // Uses scattering model to compute MFP; subclasses may override
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4CMPVScatteringRate* rateModel;	// Returns scattering rate in hertz
  G4CMPVPhysicsModel* physicsModel;	// Implements PostStepDoIt() action

  // hide assignment operators as private 
  G4CMPVProcess(G4CMPVProcess&);
  G4CMPVProcess& operator=(const G4CMPVProcess& right);
};

#endif	/* G4CMPVProcess_h */
