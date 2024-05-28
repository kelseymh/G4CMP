/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialConfigMessenger_hh
#define RISQTutorialConfigMessenger_hh 1

// $Id$
// File:  RISQTutorialConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		RISQTutorialConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class RISQTutorialConfigManager;
class G4UIcmdWithAString;
class G4UIcommand;


class RISQTutorialConfigMessenger : public G4UImessenger {
public:
  RISQTutorialConfigMessenger(RISQTutorialConfigManager* theData);
  virtual ~RISQTutorialConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RISQTutorialConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;

private:
  RISQTutorialConfigMessenger(const RISQTutorialConfigMessenger&);	// Copying is forbidden
  RISQTutorialConfigMessenger& operator=(const RISQTutorialConfigMessenger&);
};

#endif /* RISQTutorialConfigMessenger_hh */
