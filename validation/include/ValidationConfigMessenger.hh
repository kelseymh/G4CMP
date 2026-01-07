/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ValidationConfigMessenger_hh
#define ValidationConfigMessenger_hh 1

// $Id$
// File:  ValidationConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		ValidationConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class ValidationConfigManager;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;


class ValidationConfigMessenger : public G4UImessenger {
public:
  ValidationConfigMessenger(ValidationConfigManager* theData);
  virtual ~ValidationConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  ValidationConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;
  G4UIcmdWithAnInteger* geometryCmd;
  G4UIcmdWithAString* stepFileCmd;  

private:
  ValidationConfigMessenger(const ValidationConfigMessenger&);	// Copying is forbidden
  ValidationConfigMessenger& operator=(const ValidationConfigMessenger&);
};

#endif /* ValidationConfigMessenger_hh */
