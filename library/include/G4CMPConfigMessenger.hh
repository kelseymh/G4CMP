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
// 20170802  Add commands for separate Luke, downconversion scaing
// 20170815  Add command to set volume surface clearance
// 20170816  Remove directory and command handlers; G4UImessenger does it!
// 20170821  Add command to select Edelweiss IV scattering model
// 20170823  Move geometry-specific commands to examples
// 20170830  Add command to set downsampling energy scale
// 20170830  Add command to set flag for producing e/h "cloud"
// 20190711  Add command to select non-ionizing energy loss function
// 20191014  Drop command for anharmonic decay sampling.
// 20200211  Add command to report version from .g4cmp-version
// 20200411  G4CMP-195: Add commands to set charge trapping MFPs

#include "G4UImessenger.hh"

class G4CMPConfigManager;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcommand;


class G4CMPConfigMessenger : public G4UImessenger {
public:
  G4CMPConfigMessenger(G4CMPConfigManager* theData);
  virtual ~G4CMPConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  G4CMPConfigManager* theManager;

  G4UIcmdWithoutParameter* versionCmd;
  G4UIcmdWithAnInteger* verboseCmd;
  G4UIcmdWithAnInteger* ehBounceCmd;
  G4UIcmdWithAnInteger* pBounceCmd;
  G4UIcmdWithADoubleAndUnit* clearCmd;
  G4UIcmdWithADoubleAndUnit* minEPhononCmd;
  G4UIcmdWithADoubleAndUnit* minEChargeCmd;
  G4UIcmdWithADoubleAndUnit* sampleECmd;
  G4UIcmdWithADoubleAndUnit* trapEMFPCmd;
  G4UIcmdWithADoubleAndUnit* trapHMFPCmd;
  G4UIcmdWithADouble* minstepCmd;
  G4UIcmdWithADouble* makePhononCmd;
  G4UIcmdWithADouble* makeChargeCmd;
  G4UIcmdWithADouble* lukePhononCmd;
  G4UIcmdWithAString* dirCmd;
  G4UIcmdWithAString* ivRateModelCmd;
  G4UIcmdWithAString* nielPartitionCmd;
  G4UIcmdWithABool*   kvmapCmd;
  G4UIcmdWithABool*   fanoStatsCmd;
  G4UIcmdWithABool*   ehCloudCmd;

private:
  G4CMPConfigMessenger(const G4CMPConfigMessenger&);	// Copying is forbidden
  G4CMPConfigMessenger& operator=(const G4CMPConfigMessenger&);
};

#endif /* G4CMPConfigMessenger_hh */
