#ifndef G4CMPConfigMessenger_hh
#define G4CMPConfigMessenger_hh 1

// $Id$
// File:  G4CMPConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		G4CMPConfigManager.
//
// 20140904  Michael Kelsey
// 20141029  Add command to set output e/h positions file

#include "G4UImessenger.hh"

class G4CMPConfigManager;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;


class G4CMPConfigMessenger : public G4UImessenger {
public:
  G4CMPConfigMessenger(G4CMPConfigManager* theData);
  virtual ~G4CMPConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

protected:
  // Create or access directory path common to all commands
  void CreateDirectory(const char* path, const char* desc);

  // Create G4UIcommand (arbitrary subclass) within current command path
  template <class T>
  T* CreateCommand(const G4String& commandName, const G4String& description);

private:
  G4CMPConfigManager* theManager;

  G4bool localCmdDir;		// Flag if directory was created or found
  G4UIdirectory* cmdDir;

  G4UIcmdWithAnInteger* verboseCmd;
  G4UIcmdWithADoubleAndUnit* voltageCmd;
  G4UIcmdWithAString* fileCmd;
  G4UIcmdWithAString* dirCmd;
  G4UIcmdWithAString* hitsCmd;
};

#include "G4CMPConfigMessenger.icc"

#endif /* G4CMPConfigMessenger_hh */
