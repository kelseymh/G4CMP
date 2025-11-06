/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVQPProcess.hh
/// \brief Definition of the G4CMPVQPProcess base class
/// 
/// This class is a base QP class for all of the QP processes that
/// inherit from the G4VDiscreteProcess. We note that the main
/// QP diffusion process does NOT inherit from this.
//
// $Id$
//
// 20250922  G4CMP-219 -- First addition to this history (done at time
//                        of merge to develop)

#ifndef G4CMPVQPProcess_h
#define G4CMPVQPProcess_h 1

#include "G4CMPVProcess.hh"
#include "G4CMPProcessSubType.hh"

class G4ParticleDefinition;


class G4CMPVQPProcess : public G4CMPVProcess {
public:
  G4CMPVQPProcess(const G4String& processName, G4CMPProcessSubType stype)
    : G4CMPVProcess(processName, stype) {;}
  virtual ~G4CMPVQPProcess() {;}
  
  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

protected:
  virtual G4bool RandomizeFinalStateMomentumDirectionInXY();
  
private:
  // hide assignment operators as private 
  G4CMPVQPProcess(G4CMPVQPProcess&);
  G4CMPVQPProcess& operator=(const G4CMPVQPProcess& right);
};

#endif	/* G4CMPVQPProcess_h */
