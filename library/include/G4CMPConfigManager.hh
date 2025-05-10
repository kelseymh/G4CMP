/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
 \***********************************************************************/

#ifndef G4CMPConfigManager_hh
#define G4CMPConfigManager_hh 1

// $Id$
// File:  G4CMPConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		applications at runtime.  Looks for environment variables
//		at initialization to set default values; active values may
//		be changed via macro commands (see G4CMPConfigMessenger).
//
// 20140904  Michael Kelsey
// 20141231  Add parameter to set scale (relative to l0) for minimum steps
// 20150106  Move Luke phonon generating flag here, out of processes
// 20150122  Add parameter to rescale voltage in Epot field files
// 20150603  Add parameter to limit reflections in DriftBoundaryProcess
// 20160624  Add flag to use or ignore phonon KV lookup tables
// 20160830  Add parameter to scale production of e/h pairs, like phonons
// 20160901  Add parameters to set minimum energy for phonons, charges
// 20170525  Block 'rule of five' copy/move semantics, as singleton
// 20170802  Add separate scaling factors for Luke and downconversion
// 20170815  Add parameter for required clearance from volume surfaces
// 20170821  Add parameter to select Edelweiss IV scattering model
// 20170823  Remove geometry-specific parameters; implement in examples
// 20170830  Add downsampling energy scale parameter
// 20170830  Add flag to create e/h pairs in "cloud" surround location
// 20180801  Change IVEdelweiss flag to string IVRateModel.
// 20190711  G4CMP-158: Add functions to select NIEL yield functions
// 20191014  G4CMP-179: Drop sampling of anharmonic decay (downconversion)
// 20200211  G4CMP-191: Add version identification from .g4cmp-version
// 20200331  G4CMP-195: Add chage trapping MFP
// 20200331  G4CMP-196: Add impact ionization mean free path
// 20200426  G4CMP-196: Change "impact ionization" to "trap ionization"
// 20200501  G4CMP-196: Change trap-ionization MFP names, "eTrap" -> "DTrap",
//		"hTrap" -> "ATrap".
// 20200504  G4CMP-195: Reduce length of charge-trapping parameter names
// 20200530  G4CMP-202: Provide separate master and worker instances
// 20200614  G4CMP-211: Add functionality to print settings
// 20210303  G4CMP-243: Add parameter to set step length for merging hits
// 20210910  G4CMP-272: Add parameter to set number of downsampled Luke phonons
// 20220921  G4CMP-319: Add temperature setting for use with QP sensors.
// 20221117  G4CMP-343: Add option flag to preserve all internal phonons.
// 20240506  G4CMP-371: Add flag to keep or discard below-minimum track energy.
// 20241224  G4CMP-419: Add parameter to set LukeScattering debug file
// 20250209  G4CMP-457: Add short names for Lindhard empirical ionization model.
// 20250325  G4CMP-463:  Add parameter for phonon surface step size & limit.

#include "globals.hh"
#include <iosfwd>

class G4CMPConfigMessenger;
class G4VNIELPartition;


class G4CMPConfigManager {
public:
  // Thread-specific instance for use with Get/Set functions
  static G4CMPConfigManager* Instance();

  ~G4CMPConfigManager();	// Must be public for end-of-job cleanup

  // Access G4CMP as-built version number
  static const G4String& Version()      { return Instance()->version; }

  // Access G4CMP's physics ID for aux. track information
  // FIXME: This maybe should go in G4CMPVProcess when it exists.
  static G4int GetPhysicsModelID()      { return Instance()->fPhysicsModelID; }
  
  // Access current values
  static G4int GetVerboseLevel()         { return Instance()->verbose; }
  static G4int GetMaxChargeBounces()	 { return Instance()->ehBounces; }
  static G4int GetMaxPhononBounces()	 { return Instance()->pBounces; }
  static G4int GetMaxLukePhonons()       { return Instance()->maxLukePhonons; }
  static G4int GetPhononSurfStepLimit()  { return Instance()->pSurfStepLimit; }
  static G4bool UseKVSolver()            { return Instance()->useKVsolver; }
  static G4bool FanoStatisticsEnabled()  { return Instance()->fanoEnabled; }
  static G4bool KeepKaplanPhonons()      { return Instance()->kaplanKeepPh; }
  static G4bool CreateChargeCloud()      { return Instance()->chargeCloud; }
  static G4bool RecordMinETracks()       { return Instance()->recordMinE; }
  static G4double GetSurfaceClearance()  { return Instance()->clearance; }
  static G4double GetMinStepScale()      { return Instance()->stepScale; }
  static G4double GetMinPhononEnergy()   { return Instance()->EminPhonons; }
  static G4double GetMinChargeEnergy()   { return Instance()->EminCharges; }
  static G4double GetSamplingEnergy()    { return Instance()->sampleEnergy; }
  static G4double GetGenPhonons()        { return Instance()->genPhonons; }
  static G4double GetGenCharges()        { return Instance()->genCharges; }
  static G4double GetLukeSampling()      { return Instance()->lukeSample; }
  static G4double GetComboStepLength()   { return Instance()->combineSteps; }
  static G4double GetETrappingMFP()      { return Instance()->eTrapMFP; }
  static G4double GetHTrappingMFP()      { return Instance()->hTrapMFP; }
  static G4double GetEDTrapIonMFP()      { return Instance()->eDTrapIonMFP; }
  static G4double GetEATrapIonMFP()      { return Instance()->eATrapIonMFP; }
  static G4double GetHDTrapIonMFP()      { return Instance()->hDTrapIonMFP; }
  static G4double GetHATrapIonMFP()      { return Instance()->hATrapIonMFP; }
  static G4double GetTemperature()       { return Instance()->temperature; }
  static G4double GetPhononSurfStepSize()  { return Instance()->pSurfStepSize; }
  static G4double GetEmpklow()      { return Instance()->Empklow; }
  static G4double GetEmpkhigh()     { return Instance()->Empkhigh; }
  static G4double GetEmpElow()      { return Instance()->EmpElow; }
  static G4double GetEmpEhigh()     { return Instance()->EmpEhigh; }
  static G4double GetEmpkFixed()    { return Instance()->EmpkFixed; }
  static G4bool GetEmpEDepK()  { return Instance()->EmpEDepK; }


  static const G4String& GetLatticeDir() { return Instance()->LatticeDir; }
  static const G4String& GetIVRateModel() { return Instance()->IVRateModel; }
  static const G4String& GetLukeDebugFile() { return Instance()->lukeFilename; }

  static const G4VNIELPartition* GetNIELPartition() { return Instance()->nielPartition; }

  // Change values (e.g., via Messenger) -- pass strings by value for toLower()
  static void SetVerboseLevel(G4int value) { Instance()->verbose = value; }
  static void SetMaxChargeBounces(G4int value) { Instance()->ehBounces = value; }
  static void SetMaxPhononBounces(G4int value) { Instance()->pBounces = value; }
  static void SetPhononSurfStepSize(G4double value) { Instance()->pSurfStepSize = value; }
  static void SetPhononSurfStepLimit(G4int value) { Instance()->pSurfStepLimit = value; }
  static void SetMaxLukePhonons(G4int value) { Instance()->maxLukePhonons = value; }
  static void SetSurfaceClearance(G4double value) { Instance()->clearance = value; }
  static void SetMinStepScale(G4double value) { Instance()->stepScale = value; }
  static void SetMinPhononEnergy(G4double value) { Instance()->EminPhonons = value; }
  static void SetMinChargeEnergy(G4double value) { Instance()->EminCharges = value; }
  static void SetSamplingEnergy(G4double value) { Instance()->sampleEnergy = value; }
  static void SetGenPhonons(G4double value) { Instance()->genPhonons = value; }
  static void SetGenCharges(G4double value) { Instance()->genCharges = value; }
  static void SetLukeSampling(G4double value) { Instance()->lukeSample = value; }
  static void SetComboStepLength(G4double value) { Instance()->combineSteps = value; }
  static void RecordMinETracks(G4bool value) { Instance()->recordMinE = value; }
  static void UseKVSolver(G4bool value) { Instance()->useKVsolver = value; }
  static void EnableFanoStatistics(G4bool value) { Instance()->fanoEnabled = value; }
  static void KeepKaplanPhonons(G4bool value) { Instance()->kaplanKeepPh = value; }
  static void SetIVRateModel(G4String value) { Instance()->IVRateModel = value; }
  static void CreateChargeCloud(G4bool value) { Instance()->chargeCloud = value; }

  static void SetETrappingMFP(G4double value) { Instance()->eTrapMFP = value; }
  static void SetHTrappingMFP(G4double value) { Instance()->hTrapMFP = value; }
  static void SetEDTrapIonMFP(G4double value) { Instance()->eDTrapIonMFP = value; }
  static void SetEATrapIonMFP(G4double value) { Instance()->eATrapIonMFP = value; }
  static void SetHDTrapIonMFP(G4double value) { Instance()->hDTrapIonMFP = value; }
  static void SetHATrapIonMFP(G4double value) { Instance()->hATrapIonMFP = value; }
  static void SetTemperature(G4double value)  { Instance()->temperature = value; }

  static void SetLukeDebugFile(const G4String& value) { Instance()->lukeFilename = value; }

  static void SetNIELPartition(const G4String& value) { Instance()->setNIEL(value); }
  static void SetNIELPartition(G4VNIELPartition* niel) { Instance()->setNIEL(niel); }

  // Empirical Lindhard settings 
  static void SetEmpklow(G4double value) { Instance()->Empklow = value; }
  static void SetEmpkhigh(G4double value) { Instance()->Empkhigh = value; }
  static void SetEmpElow(G4double value) { Instance()->EmpElow = value; }
  static void SetEmpEhigh(G4double value) { Instance()->EmpEhigh = value; }
  static void SetEmpkFixed(G4double value) { Instance()->EmpkFixed = value; }
  static void SetEmpEDepK(G4bool value) { Instance()->EmpEDepK = value; }

  // These settings require the geometry to be rebuilt
  static void SetLatticeDir(const G4String& dir)
  { Instance()->LatticeDir=dir; UpdateGeometry(); }
  
  static void UpdateGeometry();

  // Print out all configuration settings
  static void Print(std::ostream& os) { Instance()->printConfig(os); }
  void printConfig(std::ostream& os) const;

private:
  G4CMPConfigManager();		// Singleton: only constructed in master thread
  G4CMPConfigManager(const G4CMPConfigManager&);	// To clone from master

  // Suppress moving and assignment operations; cloning is done above
  G4CMPConfigManager(G4CMPConfigManager&&) = delete;
  G4CMPConfigManager& operator=(const G4CMPConfigManager&) = delete;
  G4CMPConfigManager& operator=(G4CMPConfigManager&&) = delete;

  // Constructor will call function to read .g4cmp-version file
  void setVersion();

  // Constructor will call by-string function to map name to class
  void setNIEL(G4String value);
  void setNIEL(G4VNIELPartition* niel);

private:
  G4int verbose;	 // Global verbosity (all processes, lattices)
  G4int fPhysicsModelID; // ID key to get aux. track info.
  G4int ehBounces;	// Maximum e/h reflections ($G4CMP_EH_BOUNCES)
  G4int pBounces;	// Maximum phonon reflections ($G4CMP_PHON_BOUNCES)
  G4int maxLukePhonons; // Approx. Luke phonon limit ($G4MP_MAX_LUKE)
  G4int pSurfStepLimit;  // Phonon surface displacement step limit ($G4CMP_PHON_SURFLIMIT).
  G4String version;	 // Version name string extracted from .g4cmp-version
  G4String LatticeDir;	 // Lattice data directory ($G4LATTICEDATA)
  G4String IVRateModel;	 // Model for IV rate ($G4CMP_IV_RATE_MODEL)
  G4String lukeFilename; // Filename for LukeScattering debugging output
  G4double eTrapMFP;	 // Mean free path for electron trapping
  G4double hTrapMFP;	 // Mean free path for hole trapping
  G4double eDTrapIonMFP; // Mean free path for e- on e-trap ionization ($G4CMP_EETRAPION_MFP)
  G4double eATrapIonMFP; // Mean free path for e- on h-trap ionization ($G4CMP_EHTRAPION_MFP)
  G4double hDTrapIonMFP; // Mean free path for h+ on e-trap ionization ($G4CMP_HETRAPION_MFP)
  G4double hATrapIonMFP; // Mean free path for h+ on h-trap ionization ($G4CMP_HHTRAPION_MFP)
  G4double temperature;  // Temperature of device, substrate, sensors, etc.
  G4double clearance;	 // Minimum distance of tracks from boundaries ($G4CMP_CLEARANCE)
  G4double stepScale;	 // Fraction of l0 for steps ($G4CMP_MIN_STEP)
  G4double sampleEnergy; // Energy above which to do sampling ($G4CMP_SAMPLE_ENERGY)
  G4double genPhonons;	 // Rate to create primary phonons ($G4CMP_MAKE_PHONONS)
  G4double genCharges;	 // Rate to create primary e/h pairs ($G4CMP_MAKE_CHARGES)
  G4double lukeSample;   // Rate to create Luke phonons ($G4CMP_LUKE_SAMPLE)
  G4double combineSteps; // Maximum length to merge track steps ($G4CMP_COMBINE_STEPLEN)
  G4double EminPhonons;	 // Minimum energy to track phonons ($G4CMP_EMIN_PHONONS)
  G4double EminCharges;	 // Minimum energy to track e/h ($G4CMP_EMIN_CHARGES)
  G4double pSurfStepSize;  // Phonon surface displacement step size ($G4CMP_PHON_SURFSTEP).
  G4bool useKVsolver;	 // Use K-Vg eigensolver ($G4CMP_USE_KVSOLVER)
  G4bool fanoEnabled;	 // Apply Fano statistics to ionization energy deposits ($G4CMP_FANO_ENABLED)
  G4bool kaplanKeepPh;   // Emit or iterate over all phonons in KaplanQP ($G4CMP_KAPLAN_KEEP)
  G4bool chargeCloud;    // Produce e/h pairs around position ($G4CMP_CHARGE_CLOUD) 
  G4bool recordMinE;     // Store below-minimum track energy as NIEL when killed
  G4VNIELPartition* nielPartition; // Function class to compute non-ionizing ($G4CMP_NIEL_FUNCTION)
  // Empirical Lindhard Model Parameters
    // Model fit parameters
  G4double Empklow;  
  G4double Empkhigh; 
    // Model validity energy range
  G4double EmpElow;  
  G4double EmpEhigh;
    // Flag to use Empirical Lindhard with energy-dependent k
  G4bool EmpEDepK; 
    // If k is not energy dependent, provide/use kFixed
  G4double EmpkFixed; 
  //
  G4CMPConfigMessenger* messenger;	// User interface (UI) commands
};

// Report configuration parameter values
inline 
std::ostream& operator<<(std::ostream& os, const G4CMPConfigManager& config) {
  config.printConfig(os);
  return os;
}

#endif	/* G4CMPConfigManager_hh */
