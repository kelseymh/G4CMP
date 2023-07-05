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
// 20190711  Add command to select non-ionizing energy loss function
// 20191014  Drop command for anharmonic decay sampling.
// 20200211  Add command to report version from .g4cmp-version
// 20200411  G4CMP-195: Add commands to set charge trapping MFPs
// 20200411  G4CMP-196: Add commands to set impact ionization MFPs
// 20200426  G4CMP-196: Change "impact ionization" to "trap ionization"
// 20200501  G4CMP-196: Change trap-ionization MFP names, "eTrap" -> "DTrap",
//		"hTrap" -> "ATrap".
// 20200504  G4CMP-195:  Reduce length of charge-trapping parameter names
// 20200614  G4CMP-211:  Add functionality to print settings
// 20210303  G4CMP-243:  Add parameter to set step length for merging hits
// 20210910  G4CMP-272:  Add parameter for soft maximum Luke phonons per event
// 20220921  G4CMP-319:  Add temperature setting for use with QP sensors.
// 20221214  G4CMP-350:  Bug fix for new temperature setting units.

#include "G4CMPConfigMessenger.hh"
#include "G4CMPConfigManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"


// Constructor and destructor

G4CMPConfigMessenger::G4CMPConfigMessenger(G4CMPConfigManager* mgr)
  : G4UImessenger("/g4cmp/",
		  "User configuration for G4CMP phonon/charge carrier library"),
    theManager(mgr), versionCmd(0), printCmd(0), verboseCmd(0), ehBounceCmd(0),
    pBounceCmd(0), maxLukeCmd(0), clearCmd(0), minEPhononCmd(0),
    minEChargeCmd(0), sampleECmd(0), comboStepCmd(0), trapEMFPCmd(0),
    trapHMFPCmd(0), eDTrapIonMFPCmd(0), eATrapIonMFPCmd(0),
    hDTrapIonMFPCmd(0), hATrapIonMFPCmd(0), tempCmd(0), minstepCmd(0),
    makePhononCmd(0), makeChargeCmd(0), lukePhononCmd(0), dirCmd(0),
    ivRateModelCmd(0), nielPartitionCmd(0), kvmapCmd(0), fanoStatsCmd(0),
    ehCloudCmd(0) {
  verboseCmd = CreateCommand<G4UIcmdWithAnInteger>("verbose",
					   "Enable diagnostic messages");

  printCmd = CreateCommand<G4UIcmdWithoutParameter>("printConfig",
				    "Report G4CMP configuration settings");

  versionCmd = CreateCommand<G4UIcmdWithoutParameter>("version",
					    "Report G4CMP version string");

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

  maxLukeCmd = CreateCommand<G4UIcmdWithAnInteger>("maxLukePhonons",
		   "Set 'maximum' number of Luke phonons produced per event");
  maxLukeCmd->SetGuidance("This is a soft maximum, estimated from the bias");
  maxLukeCmd->SetGuidance("voltage of the device and the downsampling scale");

  minEPhononCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("minEPhonons",
          "Minimum energy for creating or tracking phonons");
  minEPhononCmd->SetUnitCategory("Energy");

  minEChargeCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("minECharges",
          "Minimum energy for creating or tracking charge carriers");
  minEChargeCmd->SetUnitCategory("Energy");

  comboStepCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("combiningStepLength",
	  "Maximum track step-length to merge energy deposit for partitioning");
  comboStepCmd->SetUnitCategory("Length");

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

  trapEMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("eTrappingMFP",
	   "Mean free path for trapping of electrons by D-type impurities");
  trapEMFPCmd->SetUnitCategory("Length");

  trapHMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("hTrappingMFP",
	   "Mean free path for trapping of holes by A-type impurities");
  trapHMFPCmd->SetUnitCategory("Length");

  eDTrapIonMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("eDTrapIonizationMFP",
	   "Mean free path for e-trap ionization by electrons");
  eDTrapIonMFPCmd->SetUnitCategory("Length");

  eATrapIonMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("eATrapIonizationMFP",
	   "Mean free path for h-trap ionization by electrons");
  eATrapIonMFPCmd->SetUnitCategory("Length");

  hDTrapIonMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("hDTrapIonizationMFP",
	   "Mean free path for e-trap ionization by holes");
  hDTrapIonMFPCmd->SetUnitCategory("Length");

  hATrapIonMFPCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("hATrapIonizationMFP",
	   "Mean free path for h-trap ionization by holes");
  hATrapIonMFPCmd->SetUnitCategory("Length");

  tempCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("temperature",
	   "Temperature to be used for device, substrate, sensors, etc.");
  tempCmd->SetUnitCategory("Temperature");

  nielPartitionCmd = CreateCommand<G4UIcmdWithAString>("NIELPartition",
	       "Select calculation for non-ionizing energy loss (NIEL)");
  nielPartitionCmd->SetCandidates("Lindhard lindhard Lin lin LewinSmith lewinsmith Lewin lewin Lew Lew");

  ehCloudCmd = CreateCommand<G4UIcmdWithABool>("createChargeCloud",
       "Produce e/h pairs in cloud surrounding energy deposit position");
  ehCloudCmd->SetDefaultValue(true);
}


G4CMPConfigMessenger::~G4CMPConfigMessenger() {
  delete printCmd; printCmd=0;
  delete verboseCmd; verboseCmd=0;
  delete versionCmd; versionCmd=0;
  delete ehBounceCmd; ehBounceCmd=0;
  delete pBounceCmd; pBounceCmd=0;
  delete maxLukeCmd; maxLukeCmd=0;
  delete clearCmd; clearCmd=0;
  delete minEPhononCmd; minEPhononCmd=0;
  delete minEChargeCmd; minEChargeCmd=0;
  delete sampleECmd; sampleECmd=0;
  delete comboStepCmd; comboStepCmd=0;
  delete trapEMFPCmd; trapEMFPCmd=0;
  delete trapHMFPCmd; trapHMFPCmd=0;
  delete eDTrapIonMFPCmd; eDTrapIonMFPCmd=0;
  delete eATrapIonMFPCmd; eATrapIonMFPCmd=0;
  delete hDTrapIonMFPCmd; hDTrapIonMFPCmd=0;
  delete hATrapIonMFPCmd; hATrapIonMFPCmd=0;
  delete tempCmd; tempCmd=0;
  delete minstepCmd; minstepCmd=0;
  delete makePhononCmd; makePhononCmd=0;
  delete makeChargeCmd; makeChargeCmd=0;
  delete lukePhononCmd; lukePhononCmd=0;
  delete dirCmd; dirCmd=0;
  delete kvmapCmd; kvmapCmd=0;
  delete fanoStatsCmd; fanoStatsCmd=0;
  delete ehCloudCmd; ehCloudCmd=0;
  delete ivRateModelCmd; ivRateModelCmd=0;
  delete nielPartitionCmd; nielPartitionCmd=0;
}


// Parse user input and add to configuration

void G4CMPConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == verboseCmd) theManager->SetVerboseLevel(StoI(value));
  if (cmd == minstepCmd) theManager->SetMinStepScale(StoD(value));
  if (cmd == makePhononCmd) theManager->SetGenPhonons(StoD(value));
  if (cmd == makeChargeCmd) theManager->SetGenCharges(StoD(value));
  if (cmd == lukePhononCmd) theManager->SetLukeSampling(StoD(value));
  if (cmd == maxLukeCmd) theManager->SetMaxLukePhonons(StoI(value));
  if (cmd == ehBounceCmd) theManager->SetMaxChargeBounces(StoI(value));
  if (cmd == pBounceCmd) theManager->SetMaxPhononBounces(StoI(value));
  if (cmd == dirCmd) theManager->SetLatticeDir(value);

  if (cmd == clearCmd)
    theManager->SetSurfaceClearance(clearCmd->GetNewDoubleValue(value));

  if (cmd == minEPhononCmd)
    theManager->SetMinPhononEnergy(minEPhononCmd->GetNewDoubleValue(value));

  if (cmd == minEChargeCmd)
    theManager->SetMinChargeEnergy(minEChargeCmd->GetNewDoubleValue(value));

  // TEMPORARY: If sampling energy is set and Luke=1., set Luke=-1.
  if (cmd == sampleECmd) {
    theManager->SetSamplingEnergy(sampleECmd->GetNewDoubleValue(value));
    if (theManager->GetLukeSampling() == 1.) theManager->SetLukeSampling(-1.);
  }

  if (cmd == comboStepCmd)
    theManager->SetComboStepLength(comboStepCmd->GetNewDoubleValue(value));

  if (cmd == trapEMFPCmd)
    theManager->SetETrappingMFP(trapEMFPCmd->GetNewDoubleValue(value));

  if (cmd == trapHMFPCmd)
    theManager->SetHTrappingMFP(trapHMFPCmd->GetNewDoubleValue(value));

  if (cmd == eDTrapIonMFPCmd)
    theManager->SetEDTrapIonMFP(eDTrapIonMFPCmd->GetNewDoubleValue(value));

  if (cmd == eATrapIonMFPCmd)
    theManager->SetEATrapIonMFP(eATrapIonMFPCmd->GetNewDoubleValue(value));

  if (cmd == hDTrapIonMFPCmd)
    theManager->SetHDTrapIonMFP(hDTrapIonMFPCmd->GetNewDoubleValue(value));

  if (cmd == hATrapIonMFPCmd)
    theManager->SetHATrapIonMFP(hATrapIonMFPCmd->GetNewDoubleValue(value));

  if (cmd == tempCmd)
    theManager->SetTemperature(tempCmd->GetNewDoubleValue(value));

  if (cmd == kvmapCmd) theManager->UseKVSolver(StoB(value));
  if (cmd == fanoStatsCmd) theManager->EnableFanoStatistics(StoB(value));
  if (cmd == ivRateModelCmd) theManager->SetIVRateModel(value);
  if (cmd == nielPartitionCmd) theManager->SetNIELPartition(value);
  if (cmd == ehCloudCmd) theManager->CreateChargeCloud(StoB(value));

  if (cmd == versionCmd)
    G4cout << "G4CMP version: " << theManager->Version() << G4endl;

  if (cmd == printCmd) G4cout << *theManager << G4endl;
}
