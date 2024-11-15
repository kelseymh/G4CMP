/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20241024 Israel Hernandez -- IIT, QSC and Fermilab


#ifndef Caustic_PhononConfigMessenger_hh
#define Caustic_PhononConfigMessenger_hh 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithADouble.hh"

class Caustic_PhononConfigManager;
class G4UIcmdWithAString;
class G4UIcommand;


class Caustic_PhononConfigMessenger : public G4UImessenger {
public:
  Caustic_PhononConfigMessenger(Caustic_PhononConfigManager* theData);
  virtual ~Caustic_PhononConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  Caustic_PhononConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;


private:
  Caustic_PhononConfigMessenger(const Caustic_PhononConfigMessenger&);	// Copying is forbidden
  Caustic_PhononConfigMessenger& operator=(const Caustic_PhononConfigMessenger&);
};

#endif
