/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVPhysicsModel.hh
/// \brief Abstract base class to define a physics implementation ("model")
///        associated with a process.  This will allow selection of
///        different models in the user's physics list, G4CMPPhysics, or
///        via UI commands.  Note that the constructor requires a G4CMP
///	   process be indicated as the owner.
//
// $Id$
//
// 20250430  New abstract base class, used by G4CMPVProcess.

#ifndef G4CMPVPhysicsModel_hh
#define G4CMPVPhysicsModel_hh 1

#include "G4CMPProcessUtils.hh"
#include "G4ParticleChange.hh"
#include "G4String.hh"
#include "G4Types.hh"

class G4CMPVProcess;
class G4Step;
class G4Track;


class G4CMPVPhysicsModel : public G4CMPProcessUtils {
public:
  G4CMPVPhysicsModel(const G4String& name, G4CMPVProcess* process=0)
    : G4CMPProcessUtils(), modelName(name), theProcess(process),
      verboseLevel(process?process->GetVerboseLevel():0) {;}

  virtual ~G4CMPVPhysicsModel() {;}

  // Register parent process after creation, if necessary
  void SetProcess(G4CMPVProcess* process) { theProcess = process; }
  
  // Change verbosity level at run time
  virtual void SetVerboseLevel(G4int vb=0) { verboseLevel=vb; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  // Subclasses must implement this; there is no default behaviour
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) = 0;

protected:			// Subclasses can access these for convenience
  G4String modelName;		// Name string useful for debugging messages
  G4CMPVProcess* theProcess;	// Parent process, to access e.g. rate model
  G4int verboseLevel;		// For diagnostic messages

  G4ParticleChange aParticleChange;	// Returned to process by pointer
};

#endif	/* G4CMPVPhysicsModel_hh */
