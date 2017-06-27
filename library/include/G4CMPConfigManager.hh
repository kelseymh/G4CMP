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

#include "globals.hh"
#include "G4RunManager.hh"

class G4CMPConfigMessenger;


class G4CMPConfigManager {
public:
  ~G4CMPConfigManager();	// Must be public for end-of-job cleanup

  // Access G4CMP's physics ID for aux. track information
  // FIXME: This maybe should go in G4CMPVProcess when it exists.
  static G4int GetPhysicsModelID()      { return Instance()->fPhysicsModelID; }

  // Access current values
  static G4int GetVerboseLevel()         { return Instance()->verbose; }
  static G4int GetMaxChargeBounces()	 { return Instance()->ehBounces; }
  static G4int GetMaxPhononBounces()	 { return Instance()->pBounces; }
  static G4int GetMillerH()		 { return Instance()->millerH; }
  static G4int GetMillerK()		 { return Instance()->millerK; }
  static G4int GetMillerL()		 { return Instance()->millerL; }
  static G4double GetVoltage()           { return Instance()->voltage; }
  static G4double GetMinStepScale()      { return Instance()->stepScale; }
  static G4double GetMinPhononEnergy()   { return Instance()->EminPhonons; }
  static G4double GetMinChargeEnergy()   { return Instance()->EminCharges; }
  static G4double GetGenPhonons()        { return Instance()->genPhonons; }
  static G4double GetGenCharges()        { return Instance()->genCharges; }
  static G4double GetEPotScale()         { return Instance()->epotScale; }
  static const G4String& GetEPotFile()   { return Instance()->EPot_file; }
  static const G4String& GetLatticeDir() { return Instance()->LatticeDir; }
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static G4bool UseKVSolver()            { return Instance()->useKVsolver; }
  static G4bool FanoStatisticsEnabled()  { return Instance()->fanoEnabled; }

  static void GetMillerOrientation(G4int& h, G4int& k, G4int& l) {
    h = Instance()->millerH; k = Instance()->millerK; l = Instance()->millerL;
  }

  // Change values (e.g., via Messenger)
  static void SetVerboseLevel(G4int value) { Instance()->verbose = value; }
  static void SetMaxChargeBounces(G4int value) { Instance()->ehBounces = value; }
  static void SetMaxPhononBounces(G4int value) { Instance()->pBounces = value; }
  static void SetMinStepScale(G4double value)
    { Instance()->stepScale = value; }
  static void SetMinPhononEnergy(G4double value) { Instance()->EminPhonons = value; }
  static void SetMinChargeEnergy(G4double value) { Instance()->EminCharges = value; }
  static void SetGenPhonons(G4double value) { Instance()->genPhonons = value; }
  static void SetGenCharges(G4double value) { Instance()->genCharges = value; }
  static void UseKVSolver(G4bool value) { Instance()->useKVsolver = value; }
  static void EnableFanoStatistics(G4bool value) { Instance()->fanoEnabled = value; }

  // These settings require the geometry to be rebuilt
  static void SetVoltage(G4double value)
    { Instance()->voltage = value; UpdateGeometry(); }
  static void SetEPotScale(G4double value)
    { Instance()->epotScale = value; UpdateGeometry(); }
  static void SetEPotFile(const G4String& name)
    { Instance()->EPot_file=name; UpdateGeometry(); }
  static void SetLatticeDir(const G4String& dir)
    { Instance()->LatticeDir=dir; UpdateGeometry(); }
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  static void SetMillerOrientation(G4int h, G4int k, G4int l)
    { Instance()->millerH=h; Instance()->millerK=k, Instance()->millerL=l;
      UpdateGeometry();
    }

  static void UpdateGeometry()
    { G4RunManager::GetRunManager()->ReinitializeGeometry(true); }

private:
  G4CMPConfigManager();		// Singleton: only constructed on request
  G4CMPConfigManager(const G4CMPConfigManager&) = delete;
  G4CMPConfigManager(G4CMPConfigManager&&) = delete;
  G4CMPConfigManager& operator=(const G4CMPConfigManager&) = delete;
  G4CMPConfigManager& operator=(G4CMPConfigManager&&) = delete;

  static G4CMPConfigManager* Instance();   // Only needed by static accessors
  static G4CMPConfigManager* theInstance;

private:
  G4double voltage;	// Uniform field voltage ($G4CMP_VOLTAGE)
  G4double stepScale;	// Fraction of l0 for steps ($G4CMP_MIN_STEP)
  G4double genPhonons;	// Rate to create phonons ($G4CMP_MAKE_PHONONS)
  G4double genCharges;	// Rate to create e/h pairs ($G4CMP_MAKE_CHARGES)
  G4double EminPhonons;	// Minimum energy to track phonons ($G4CMP_EMIN_PHONONS)
  G4double EminCharges;	// Minimum energy to track e/h ($G4CMP_EMIN_CHARGES)
  G4double epotScale;	// Scale factor for EPot ($G4CMP_EPOT_SCALE)
  G4int verbose;	// Global verbosity (all processes, lattices)
  G4int fPhysicsModelID; // ID key to get aux. track info.
  G4int ehBounces;	// Maximum e/h reflections ($G4CMP_EH_BOUNCES)
  G4int pBounces;	// Maximum phonon reflections ($G4CMP_PHON_BOUNCES)
  G4int millerH;	// Lattice orientation ($G4CMP_MILLER_H,_K,_L)
  G4int millerK;
  G4int millerL;
  G4bool useKVsolver;	// Use K-Vg eigensolver ($G4CMP_USE_KVSOLVER)
  G4bool fanoEnabled;	// Apply Fano statistics to ionization energy deposits
                        // ($G4CMP_FANO_ENABLED)
  G4String EPot_file;	// Name of E-field file ($G4CMP_EPOT_FILE)
  G4String LatticeDir;	// Lattice data directory ($G4LATTICEDATA)
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)

  G4CMPConfigMessenger* messenger;
};

#endif	/* G4CMPConfigManager_hh */
