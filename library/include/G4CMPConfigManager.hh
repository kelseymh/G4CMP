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
// 20170823  Remove geometry-specific parameters; implement in examples
// 20170830  Add downsampling energy scale parameter

#include "globals.hh"

class G4CMPConfigMessenger;


class G4CMPConfigManager {
public:
  static G4CMPConfigManager* Instance();
  ~G4CMPConfigManager();	// Must be public for end-of-job cleanup

  // Access G4CMP's physics ID for aux. track information
  // FIXME: This maybe should go in G4CMPVProcess when it exists.
  static G4int GetPhysicsModelID()      { return Instance()->fPhysicsModelID; }

  // Access current values
  static G4int GetVerboseLevel()         { return Instance()->verbose; }
  static G4int GetMaxChargeBounces()	 { return Instance()->ehBounces; }
  static G4int GetMaxPhononBounces()	 { return Instance()->pBounces; }
  static G4bool UseKVSolver()            { return Instance()->useKVsolver; }
  static G4bool FanoStatisticsEnabled()  { return Instance()->fanoEnabled; }
  static G4double GetSurfaceClearance()  { return Instance()->clearance; }
  static G4double GetMinStepScale()      { return Instance()->stepScale; }
  static G4double GetMinPhononEnergy()   { return Instance()->EminPhonons; }
  static G4double GetMinChargeEnergy()   { return Instance()->EminCharges; }
  static G4double GetSamplingEnergy()    { return Instance()->sampleEnergy; }
  static G4double GetGenPhonons()        { return Instance()->genPhonons; }
  static G4double GetGenCharges()        { return Instance()->genCharges; }
  static G4double GetLukeSampling()      { return Instance()->lukeSample; }
  static G4double GetDownconversionSampling() { return Instance()->downSample; }
  static const G4String& GetLatticeDir() { return Instance()->LatticeDir; }

  // Change values (e.g., via Messenger)
  static void SetVerboseLevel(G4int value) { Instance()->verbose = value; }
  static void SetMaxChargeBounces(G4int value) { Instance()->ehBounces = value; }
  static void SetMaxPhononBounces(G4int value) { Instance()->pBounces = value; }
  static void SetSurfaceClearance(G4double value) { Instance()->clearance = value; }
  static void SetMinStepScale(G4double value) { Instance()->stepScale = value; }
  static void SetMinPhononEnergy(G4double value) { Instance()->EminPhonons = value; }
  static void SetMinChargeEnergy(G4double value) { Instance()->EminCharges = value; }
  static void SetSamplingEnergy(G4double value) { Instance()->sampleEnergy = value; }
  static void SetGenPhonons(G4double value) { Instance()->genPhonons = value; }
  static void SetGenCharges(G4double value) { Instance()->genCharges = value; }
  static void SetLukeSampling(G4double value) { Instance()->lukeSample = value; }
  static void SetDownconversionSampling(G4double value) { Instance()->downSample = value; }
  static void UseKVSolver(G4bool value) { Instance()->useKVsolver = value; }
  static void EnableFanoStatistics(G4bool value) { Instance()->fanoEnabled = value; }

  // These settings require the geometry to be rebuilt
  static void SetLatticeDir(const G4String& dir)
    { Instance()->LatticeDir=dir; UpdateGeometry(); }

  static void UpdateGeometry();

private:
  G4CMPConfigManager();		// Singleton: only constructed on request
  G4CMPConfigManager(const G4CMPConfigManager&) = delete;
  G4CMPConfigManager(G4CMPConfigManager&&) = delete;
  G4CMPConfigManager& operator=(const G4CMPConfigManager&) = delete;
  G4CMPConfigManager& operator=(G4CMPConfigManager&&) = delete;

  static G4CMPConfigManager* theInstance;

private:
  G4int verbose;	// Global verbosity (all processes, lattices)
  G4int fPhysicsModelID; // ID key to get aux. track info.
  G4int ehBounces;	// Maximum e/h reflections ($G4CMP_EH_BOUNCES)
  G4int pBounces;	// Maximum phonon reflections ($G4CMP_PHON_BOUNCES)
  G4String LatticeDir;	// Lattice data directory ($G4LATTICEDATA)
  G4double clearance;	// Minimum distance of tracks from boundaries ($G4CMP_CLEARANCE)
  G4double stepScale;	// Fraction of l0 for steps ($G4CMP_MIN_STEP)
  G4double sampleEnergy; // Energy above which to do sampling ($G4CMP_SAMPLE_ENERGY)
  G4double genPhonons;	// Rate to create primary phonons ($G4CMP_MAKE_PHONONS)
  G4double genCharges;	// Rate to create primary e/h pairs ($G4CMP_MAKE_CHARGES)
  G4double lukeSample;  // Rate to create Luke phonons ($G4CMP_LUKE_SAMPLE)
  G4double downSample;  // Rate to apply downconversion ($G4CMP_DOWN_SAMPLE)
  G4double EminPhonons;	// Minimum energy to track phonons ($G4CMP_EMIN_PHONONS)
  G4double EminCharges;	// Minimum energy to track e/h ($G4CMP_EMIN_CHARGES)
  G4bool useKVsolver;	// Use K-Vg eigensolver ($G4CMP_USE_KVSOLVER)
  G4bool fanoEnabled;	// Apply Fano statistics to ionization energy deposits
                        // ($G4CMP_FANO_ENABLED)

  G4CMPConfigMessenger* messenger;
};

#endif	/* G4CMPConfigManager_hh */
