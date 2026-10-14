/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file  HeterostructureConfigMessenger.hh
/// \brief Definition of the HeterostructureConfigMessenger class
///   Macro command to set user configuration in
///   HeterostructureConfigManager.

//    20260815 Selby Q. Dang

#ifndef HeterostructureConfigMessenger_hh
#define HeterostructureConfigMessenger_hh 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithADouble.hh"  

class HeterostructureConfigManager;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcommand;


class HeterostructureConfigMessenger : public G4UImessenger {
public:
  HeterostructureConfigMessenger(HeterostructureConfigManager* theData);
  virtual ~HeterostructureConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  HeterostructureConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;
  G4UIcmdWithADouble* qAbsProbGeSiCmd;
  G4UIcmdWithADouble* qReflProbGeSiCmd;
  G4UIcmdWithADouble* qAbsProbSiGeCmd;
  G4UIcmdWithADouble* qReflProbSiGeCmd;

private:
  HeterostructureConfigMessenger(const HeterostructureConfigMessenger&);	// Copying is forbidden
  HeterostructureConfigMessenger& operator=(const HeterostructureConfigMessenger&);
};

#endif /* HeterostructureConfigMessenger_hh */
