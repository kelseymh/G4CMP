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
// 20160517  

#include "globals.hh"
#include "G4RunManager.hh"

class G4CMPConfigMessenger;


class G4CMPConfigManager {
public:
  ~G4CMPConfigManager();	// Must be public for end-of-job cleanup

  // Access current values
  static G4int GetVerboseLevel()         { return Instance()->verbose; }
  static G4int GetMaxChargeBounces()	 { return Instance()->ehBounces; }
  static G4int GetMaxPhononBounces()	 { return Instance()->pBounces; }
  static G4int GetMillerH()		 { return Instance()->millerH; }
  static G4int GetMillerK()		 { return Instance()->millerK; }
  static G4int GetMillerL()		 { return Instance()->millerL; }
  static G4double GetVoltage()           { return Instance()->voltage; }
  static G4double GetMinStepScale()      { return Instance()->stepScale; }
  static G4double GetGenPhonons()        { return Instance()->genPhonons; }
  static G4double GetEpotScale()         { return Instance()->epotScale; }
  static const G4String& GetEpotFile()   { return Instance()->Epot_file; }
  static const G4String& GetLatticeDir() { return Instance()->LatticeDir; }
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }

  static void GetMillerOrientation(G4int& h, G4int& k, G4int& l) {
    h = Instance()->millerH; k = Instance()->millerK; l = Instance()->millerL;
  }

  // Change values (e.g., via Messenger)
  static void SetVerboseLevel(G4int value)
    { Instance()->verbose = value; }
  static void SetMaxChargeBounces(G4int value)
    { Instance()->ehBounces = value; }
  static void SetMaxPhononBounces(G4int value)
    { Instance()->pBounces = value; }
  static void SetMillerOrientation(G4int h, G4int k, G4int l)
    { Instance()->millerH=h; Instance()->millerK=k, Instance()->millerL=l; }
  static void SetVoltage(G4double value)
    { Instance()->voltage = value; UpdateGeometry(); }
  static void SetMinStepScale(G4double value)
    { Instance()->stepScale = value; }
  static void SetGenPhonons(G4double value)
    { Instance()->genPhonons=value; }
  static void SetEpotScale(G4double value)
    { Instance()->epotScale = value; UpdateGeometry(); }
  static void SetEpotFile(const G4String& name)
    { Instance()->Epot_file=name; UpdateGeometry(); }
  static void SetLatticeDir(const G4String& dir)
    { Instance()->LatticeDir=dir; UpdateGeometry(); }
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  static void UpdateGeometry()
    { G4RunManager::GetRunManager()->ReinitializeGeometry(true); }

private:
  G4CMPConfigManager();		// Singleton: only constructed on request

  static G4CMPConfigManager* Instance();   // Only needed by static accessors
  static G4CMPConfigManager* theInstance;

private:
  G4double voltage;		// Uniform field voltage ($G4CMP_VOLTAGE)
  G4double stepScale;		// Fraction of l0 for steps ($G4CMP_MIN_STEP)
  G4double genPhonons;         // Rate to create phonons ($G4CMP_LUKE_PHONONS)
  G4double epotScale;		// Scale factor for Epot ($G4CMP_EPOT_SCALE)
  G4int verbose;		// Global verbosity (all processes, lattices)
  G4int ehBounces;		// Maximum e/h reflections ($G4CMP_EH_BOUNCES)
  G4int pBounces;		// Maximum phonon reflections ($G4CMP_PHON_BOUNCES)
  G4int millerH;		// Lattice orientation ($G4CMP_MILLER_H,_K,_L)
  G4int millerK;
  G4int millerL;
  G4String Epot_file;		// Name of E-field file ($G4CMP_EPOT_FILE)
  G4String LatticeDir;		// Lattice data directory ($G4LATTICEDATA)
  G4String Hit_file;		// Output file of e/h hits ($G4CMP_HIT_FILE)

  G4CMPConfigMessenger* messenger;
};

#endif	/* G4CMPConfigManager_hh */
