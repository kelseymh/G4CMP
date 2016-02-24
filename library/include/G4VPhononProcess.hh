/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4VPhononProcess.hh
/// \brief Definition of the G4VPhononProcess base class
//
// $Id$
//
// 20140312  Move utility functions to separate class, multiple inheritance
// 20140331  Add required subtype code to constructor

#ifndef G4VPhononProcess_h
#define G4VPhononProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"

class G4PhononTrackMap;
class G4LatticePhysical;


class G4VPhononProcess : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4VPhononProcess(const G4String& processName, G4CMPProcessSubType stype);
  virtual ~G4VPhononProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

  // Initialize wave vectors for currently active track(s)
  // NOTE:  These functions must call back to base class implementations!
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

private:
  // hide assignment operators as private 
  G4VPhononProcess(G4VPhononProcess&);
  G4VPhononProcess& operator=(const G4VPhononProcess& right);
};

#endif	/* G4VPhononProcess_h */
