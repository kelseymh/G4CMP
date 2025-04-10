/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleConfigMessenger_hh
#define QuasiparticleConfigMessenger_hh 1

// $Id$
// File:  QuasiparticleConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		QuasiparticleConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class QuasiparticleConfigManager;
class G4UIcmdWithAString;
class G4UIcommand;


class QuasiparticleConfigMessenger : public G4UImessenger {
public:
  QuasiparticleConfigMessenger(QuasiparticleConfigManager* theData);
  virtual ~QuasiparticleConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  QuasiparticleConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;

private:
  QuasiparticleConfigMessenger(const QuasiparticleConfigMessenger&);	// Copying is forbidden
  QuasiparticleConfigMessenger& operator=(const QuasiparticleConfigMessenger&);
};

#endif /* QuasiparticleConfigMessenger_hh */
