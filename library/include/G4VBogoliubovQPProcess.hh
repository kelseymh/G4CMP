/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4VBogoliubovQPProcess.hh
/// \brief Definition of the G4VBogoliubovQPProcess base class
//
// $Id$
//
// 20140312  Move utility functions to separate class, multiple inheritance
// 20140331  Add required subtype code to constructor
// 20170601  Inherit from new G4CMPVProcess, only need IsApplicable

#ifndef G4VBogoliubovQPProcess_h
#define G4VBogoliubovQPProcess_h 1

#include "G4CMPVProcess.hh"
#include "G4CMPProcessSubType.hh"

class G4ParticleDefinition;


class G4VBogoliubovQPProcess : public G4CMPVProcess {
public:
  G4VBogoliubovQPProcess(const G4String& processName, G4CMPProcessSubType stype)
    : G4CMPVProcess(processName, stype) {;}
  virtual ~G4VBogoliubovQPProcess() {;}
  
  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);
  virtual G4bool RandomizeFinalStateMomentumDirectionInXY();
  
private:
  // hide assignment operators as private 
  G4VBogoliubovQPProcess(G4VBogoliubovQPProcess&);
  G4VBogoliubovQPProcess& operator=(const G4VBogoliubovQPProcess& right);
};

#endif	/* G4VBogoliubovQPProcess_h */
