/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChargeConfigMessenger_hh
#define ChargeConfigMessenger_hh 1

// $Id$
// File:  ChargeConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		ChargeConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class ChargeConfigManager;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcommand;


class ChargeConfigMessenger : public G4UImessenger {
public:
  ChargeConfigMessenger(ChargeConfigManager* theData);
  virtual ~ChargeConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  ChargeConfigManager* theManager;
  G4UIcmdWithADoubleAndUnit* voltageCmd;
  G4UIcmdWithADouble* escaleCmd;
  G4UIcmdWithAString* fileCmd;
  G4UIcmdWithAString* hitsCmd;
  G4UIcmdWithAString* millerCmd;	// Will parse out three integers

private:
  ChargeConfigMessenger(const ChargeConfigMessenger&);	// Copying is forbidden
  ChargeConfigMessenger& operator=(const ChargeConfigMessenger&);
};

#endif /* ChargeConfigMessenger_hh */
