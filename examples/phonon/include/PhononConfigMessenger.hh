/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef PhononConfigMessenger_hh
#define PhononConfigMessenger_hh 1

// $Id$
// File:  PhononConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		PhononConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class PhononConfigManager;
class G4UIcmdWithAString;
class G4UIcommand;


class PhononConfigMessenger : public G4UImessenger {
public:
  PhononConfigMessenger(PhononConfigManager* theData);
  virtual ~PhononConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  PhononConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;

private:
  PhononConfigMessenger(const PhononConfigMessenger&);	// Copying is forbidden
  PhononConfigMessenger& operator=(const PhononConfigMessenger&);
};

#endif /* PhononConfigMessenger_hh */
