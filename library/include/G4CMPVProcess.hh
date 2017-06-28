/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVProcess.hh
/// \brief Definition of the G4CMPVProcess base class
//
// $Id$
//
// 20170601  New abstract base class for all G4CMP processes

#ifndef G4CMPVProcess_h
#define G4CMPVProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPProcessUtils.hh"


class G4CMPVProcess : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4CMPVProcess(const G4String& processName, G4CMPProcessSubType stype);
  virtual ~G4CMPVProcess() {;}

  // Initialize track/volume information (lattice, wavevector, etc.)
  // NOTE:  Subclasses must call back to these base implementations!
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

private:
  // hide assignment operators as private 
  G4CMPVProcess(G4CMPVProcess&);
  G4CMPVProcess& operator=(const G4CMPVProcess& right);
};

#endif	/* G4CMPVProcess_h */
