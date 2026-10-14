/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// /file  HeterostructureConfigMessenger.cc
/// /brief Implementation of the HeterostructureConfigMessenger class
///   Macro command defitions to set user configuration in
///		HeterostructureConfigManager.

// 20260815  Selby Q. Dang (Stanford/SLAC)

#include "HeterostructureConfigMessenger.hh"
#include "HeterostructureConfigManager.hh"
#include "G4UIcmdWithAString.hh"


// Constructor and destructor

HeterostructureConfigMessenger::HeterostructureConfigMessenger(HeterostructureConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP Heterostructure example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
			      "Set filename for output of Heterostructure hit locations");

  qAbsProbGeSiCmd = CreateCommand<G4UIcmdWithADouble> ("qAbsProbGeSi", 
            "Set the qAbsProb parameter from Ge to Si, between 0 and 1");
  qAbsProbGeSiCmd->SetParameterName("qAbsProbGeSi", false);
  qAbsProbGeSiCmd->SetRange("qAbsProbGeSi >= 0. && qAbsProbGeSi <= 1.");

  qAbsProbSiGeCmd = CreateCommand<G4UIcmdWithADouble> ("qAbsProbSiGe", 
            "Set the qAbsProb parameter from Si to Ge, between 0 and 1");
  qAbsProbSiGeCmd->SetParameterName("qAbsProbSiGe", false);
  qAbsProbSiGeCmd->SetRange("qAbsProbSiGe >= 0. && qAbsProbSiGe <= 1.");

  qReflProbGeSiCmd = CreateCommand<G4UIcmdWithADouble> ("qReflProbGeSi", 
            "Set the pReflProb parameter from Ge to Si, between 0 and 1");
  qReflProbGeSiCmd->SetParameterName("value", false);
  qReflProbGeSiCmd->SetRange("value >= 0. && value <= 1.");

  qReflProbSiGeCmd = CreateCommand<G4UIcmdWithADouble> ("qReflProbSiGe", 
            "Set the pReflProb parameter from Si to Ge, between 0 and 1");
  qReflProbSiGeCmd->SetParameterName("value", false);
  qReflProbSiGeCmd->SetRange("value >= 0. && value <= 1.");
}


HeterostructureConfigMessenger::~HeterostructureConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
  delete qAbsProbGeSiCmd; qAbsProbGeSiCmd=0;
  delete qAbsProbSiGeCmd; qAbsProbSiGeCmd=0;
  delete qReflProbGeSiCmd; qReflProbGeSiCmd=0;
  delete qReflProbSiGeCmd; qReflProbSiGeCmd=0;
}


// Parse user input and add to configuration

void HeterostructureConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);
  if (cmd == qAbsProbGeSiCmd) theManager->SetQAbsProbGeSi(StoD(value));
  if (cmd == qAbsProbSiGeCmd) theManager->SetQAbsProbSiGe(StoD(value));
  if (cmd == qReflProbGeSiCmd) theManager->SetQReflProbGeSi(StoD(value));
  if (cmd == qReflProbSiGeCmd) theManager->SetQReflProbSiGe(StoD(value));
}
