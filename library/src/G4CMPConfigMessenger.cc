/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  G4CMPConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		G4CMPConfigManager.
//
// 20140904  Michael Kelsey
// 20141029  Add command to set output e/h positions file
// 20150106  Add command to toggle generate Luke phonons
// 20150122  Add command to rescale Epot file voltage by some factor
// 20150603  Add command to limit reflections in DriftBoundaryProcess
// 20160518  Add commands for Miller orientation, phonon bounces
// 20160624  Add command to select KV lookup tables vs. calculator
// 20160830  Add command to scale production of e/h pairs, like phonons
// 20170802  Add commands for separate Luke, downconversion scaing
// 20170815  Add command to set volume surface clearance
// 20170816  Remove directory and command handlers; G4UImessenger does it!
// 20170821  Add command to select Edelweiss IV scattering model
// 20170823  Move geometry-specific commands to examples
// 20170830  Add command for downsampling energy scale parameter
// 20170830  Add command to set flag for producing e/h "cloud"

#include "G4CMPConfigMessenger.hh"
#include "G4CMPConfigManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


// Constructor and destructor

G4CMPConfigMessenger::G4CMPConfigMessenger(G4CMPConfigManager* mgr)
  : G4UImessenger("/g4cmp/",
		  "User configuration for G4CMP phonon/charge carrier library"),
    theManager(mgr), verboseCmd(0), ehBounceCmd(0), pBounceCmd(0), clearCmd(0),
    minEPhononCmd(0), minEChargeCmd(0), minstepCmd(0), makePhononCmd(0),
    makeChargeCmd(0), lukePhononCmd(0), downconvCmd(0),
    dirCmd(0), ivRateModelCmd(0), kvmapCmd(0), fanoStatsCmd(0), ehCloudCmd(0) {
  verboseCmd = CreateCommand<G4UIcmdWithAnInteger>("verbose",
					   "Enable diagnostic messages");

  dirCmd = CreateCommand<G4UIcmdWithAString>("LatticeData",
			     "Set directory for lattice configuration files");
  dirCmd->AvailableForStates(G4State_PreInit);

  clearCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("clearance",
	      "Minimum distance from volume boundaries for new tracks");
  clearCmd->SetUnitCategory("Length");

  minstepCmd = CreateCommand<G4UIcmdWithADouble>("minimumStep",
			 "Set fraction of L0 for charge carrier minimum step");

  sampleECmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("samplingEnergy",
			"Energy scale above which events/hits are downsampled");
  sampleECmd->SetGuidance("Below this energy, Geant4 energy deposits are");
  sampleECmd->SetGuidance("fully converted to charge carriers and phonons");
  sampleECmd->SetGuidance("by EnergyPartition.  Above, the conversion is");
  sampleECmd->SetGuidance("scaled by the ratio of this parameter to the G4");
  sampleECmd->SetGuidance("energy deposit.  This parameter overrides the");
  sampleECmd->SetGuidance("sampling rates 'producePhonons', 'produceCharges',");
  sampleECmd->SetGuidance("and 'sampleLuke'.");
  sampleECmd->SetUnitCategory("Energy");

  makePhononCmd = CreateCommand<G4UIcmdWithADouble>("producePhonons",
		    "Set rate of production of primary phonons");

  makeChargeCmd = CreateCommand<G4UIcmdWithADouble>("produceCharges",
		    "Set rate of production of primary charge carriers");

  lukePhononCmd = CreateCommand<G4UIcmdWithADouble>("sampleLuke",
		    "Set rate of Luke actual phonon production");

  downconvCmd = CreateCommand<G4UIcmdWithADouble>("downconvertPhonons",
		  "Set scale factor for rate of phonon downconversion process");

  minEPhononCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("minEPhonons",
          "Minimum energy for creating or tracking phonons");

  minEChargeCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("minECharges",
          "Minimum energy for creating or tracking charge carriers");

  ehBounceCmd = CreateCommand<G4UIcmdWithAnInteger>("chargeBounces",
		  "Maximum number of reflections allowed for charge carriers");

  pBounceCmd = CreateCommand<G4UIcmdWithAnInteger>("phononBounces",
		  "Maximum number of reflections allowed for phonons");

  kvmapCmd = CreateCommand<G4UIcmdWithABool>("useKVsolver",
			     "Use eigenvector solver for K-Vg conversion");
  kvmapCmd->SetParameterName("lookup",true,false);
  kvmapCmd->SetDefaultValue(true);

  fanoStatsCmd = CreateCommand<G4UIcmdWithABool>("enableFanoStatistics",
           "Modify input ionization energy according to Fano statistics.");
  fanoStatsCmd->SetDefaultValue(true);

  ivRateModelCmd = CreateCommand<G4UIcmdWithAString>("IVRateModel",
           "Set the model for IV scattering rate.");
  ivRateModelCmd->SetGuidance("IVRate	  : Scattering matrix calculation");
  ivRateModelCmd->SetGuidance("Linear	  : Gamma0 + Gamma * E^x");
  ivRateModelCmd->SetGuidance("Quadratic  : Gamma * sqrt[(E0^2 + E^2)^x]");
  ivRateModelCmd->SetCandidates("IVRate Linear Quadratic");
  ivRateModelCmd->SetDefaultValue("Quadratic");

  ehCloudCmd = CreateCommand<G4UIcmdWithABool>("createChargeCloud",
       "Produce e/h pairs in cloud surrounding energy deposit position");
  ehCloudCmd->SetDefaultValue(true);
}


G4CMPConfigMessenger::~G4CMPConfigMessenger() {
  delete verboseCmd; verboseCmd=0;
  delete ehBounceCmd; ehBounceCmd=0;
  delete pBounceCmd; pBounceCmd=0;
  delete clearCmd; clearCmd=0;
  delete minEPhononCmd; minEPhononCmd=0;
  delete minEChargeCmd; minEChargeCmd=0;
  delete sampleECmd; sampleECmd=0;
  delete minstepCmd; minstepCmd=0;
  delete makePhononCmd; makePhononCmd=0;
  delete makeChargeCmd; makeChargeCmd=0;
  delete lukePhononCmd; lukePhononCmd=0;
  delete downconvCmd; downconvCmd=0;
  delete dirCmd; dirCmd=0;
  delete kvmapCmd; kvmapCmd=0;
  delete fanoStatsCmd; fanoStatsCmd=0;
  delete ehCloudCmd; ehCloudCmd=0;
  delete ivRateModelCmd; ivRateModelCmd=0;
}


// Parse user input and add to configuration

void G4CMPConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == verboseCmd) theManager->SetVerboseLevel(StoI(value));
  if (cmd == minstepCmd) theManager->SetMinStepScale(StoD(value));
  if (cmd == makePhononCmd) theManager->SetGenPhonons(StoD(value));
  if (cmd == makeChargeCmd) theManager->SetGenCharges(StoD(value));
  if (cmd == lukePhononCmd) theManager->SetLukeSampling(StoD(value));
  if (cmd == downconvCmd) theManager->SetDownconversionSampling(StoD(value));
  if (cmd == ehBounceCmd) theManager->SetMaxChargeBounces(StoI(value));
  if (cmd == pBounceCmd) theManager->SetMaxPhononBounces(StoI(value));
  if (cmd == dirCmd) theManager->SetLatticeDir(value);

  if (cmd == clearCmd)
    theManager->SetSurfaceClearance(clearCmd->GetNewDoubleValue(value));

  if (cmd == minEPhononCmd)
    theManager->SetMinPhononEnergy(minEPhononCmd->GetNewDoubleValue(value));

  if (cmd == minEChargeCmd)
    theManager->SetMinChargeEnergy(minEChargeCmd->GetNewDoubleValue(value));

  if (cmd == sampleECmd)
    theManager->SetSamplingEnergy(sampleECmd->GetNewDoubleValue(value));

  if (cmd == kvmapCmd) theManager->UseKVSolver(StoB(value));
  if (cmd == fanoStatsCmd) theManager->EnableFanoStatistics(StoB(value));
  if (cmd == ivRateModelCmd) theManager->SetIVRateModel(value);
  if (cmd == ehCloudCmd) theManager->CreateChargeCloud(StoB(value));
}
