/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

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
// 20141231  Add command to set scale (relative to l0) for minimum steps
// 20150106  Add command to toggle generating Luke phonons
// 20150122  Add command to rescale Epot file voltage by some factor
// 20150603  Add command to limit reflections in DriftBoundaryProcess
// 20160518  Add commands for Miller orientation, phonon bounces
// 20160624  Add command to select KV lookup tables vs. calculator
// 20160830  Add command to scale production of e/h pairs, like phonons
// 20160901  Add commands to set minimum energy for phonons, charges
// 20170821  Add command to select Edelweiss IV scattering model

#include "G4UImessenger.hh"

class G4CMPConfigManager;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
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
  G4UIcmdWithAnInteger* ehBounceCmd;
  G4UIcmdWithAnInteger* pBounceCmd;
  G4UIcmdWithADoubleAndUnit* voltageCmd;
  G4UIcmdWithADoubleAndUnit* minEPhononCmd;
  G4UIcmdWithADoubleAndUnit* minEChargeCmd;
  G4UIcmdWithADouble* minstepCmd;
  G4UIcmdWithADouble* makePhononCmd;
  G4UIcmdWithADouble* makeChargeCmd;
  G4UIcmdWithADouble* escaleCmd;
  G4UIcmdWithAString* fileCmd;
  G4UIcmdWithAString* dirCmd;
  G4UIcmdWithAString* hitsCmd;
  G4UIcmdWithAString* millerCmd;	// Will parse out three integers
  G4UIcmdWithABool*   kvmapCmd;
  G4UIcmdWithABool*   fanoStatsCmd;
  G4UIcmdWithABool*   ivEdelCmd;

private:
  G4CMPConfigMessenger(const G4CMPConfigMessenger&);	// Copying is forbidden
  G4CMPConfigMessenger& operator=(const G4CMPConfigMessenger&);
};

#include "G4CMPConfigMessenger.icc"

#endif /* G4CMPConfigMessenger_hh */
