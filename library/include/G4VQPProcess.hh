/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4VQPProcess.hh
/// \brief Definition of the G4VQPProcess base class
//
// $Id$
//
// 20140312  Move utility functions to separate class, multiple inheritance
// 20140331  Add required subtype code to constructor
// 20170601  Inherit from new G4CMPVProcess, only need IsApplicable

#ifndef G4VQPProcess_h
#define G4VQPProcess_h 1

#include "G4CMPVProcess.hh"
#include "G4CMPProcessSubType.hh"

class G4ParticleDefinition;


class G4VQPProcess : public G4CMPVProcess {
public:
  G4VQPProcess(const G4String& processName, G4CMPProcessSubType stype)
    : G4CMPVProcess(processName, stype) {;}
  virtual ~G4VQPProcess() {;}
  
  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);
  virtual G4bool RandomizeFinalStateMomentumDirectionInXY();
  
private:
  // hide assignment operators as private 
  G4VQPProcess(G4VQPProcess&);
  G4VQPProcess& operator=(const G4VQPProcess& right);
};

#endif	/* G4VQPProcess_h */
