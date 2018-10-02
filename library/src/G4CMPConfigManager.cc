/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  G4CMPConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		applications at runtime.  Looks for environment variables
//		at initialization to set default values; active values may
//		be changed via macro commands (see G4CMPConfigMessenger).
//
// 20140904  Michael Kelsey
// 20141029  Force numerical voltage to correct units
// 20150603  Add parameter to limit reflections in DriftBoundaryProcess
// 20161028  Drop default filename for EPot (mesh) field
// 20170802  Add separate scaling factors for Luke and downconversion
// 20170815  Add parameter for required clearance from volume surfaces
// 20170823  Remove geometry-specific parameters; implement in examples
// 20170830  Add downsampling energy scale parameter
// 20170908  G4CMP-118:  Use Edelweiss IV rate by default

#include "G4CMPConfigManager.hh"
#include "G4CMPConfigMessenger.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

G4CMPConfigManager* G4CMPConfigManager::theInstance = 0;

G4CMPConfigManager* G4CMPConfigManager::Instance() {
  if (!theInstance) theInstance = new G4CMPConfigManager;
  return theInstance;
}

G4CMPConfigManager::G4CMPConfigManager()
  : verbose(getenv("G4CMP_DEBUG")?atoi(getenv("G4CMP_DEBUG")):0),
    ehBounces(getenv("G4CMP_EH_BOUNCES")?atoi(getenv("G4CMP_EH_BOUNCES")):1),
    pBounces(getenv("G4CMP_PHON_BOUNCES")?atoi(getenv("G4CMP_PHON_BOUNCES")):100),
    LatticeDir(getenv("G4LATTICEDATA")?getenv("G4LATTICEDATA"):"./CrystalMaps"),
    IVRateModel(getenv("G4CMP_IV_RATE_MODEL")?getenv("G4CMP_IV_RATE_MODEL"):"Quadratic"),
    clearance(getenv("G4CMP_CLEARANCE")?strtod(getenv("G4CMP_CLEARANCE"),0)*mm:1e-6*mm),
    stepScale(getenv("G4CMP_MIN_STEP")?strtod(getenv("G4CMP_MIN_STEP"),0):-1.),
    sampleEnergy(getenv("G4CMP_SAMPLE_ENERGY")?strtod(getenv("G4CMP_SAMPLE_ENERGY"),0):-1.),
    genPhonons(getenv("G4CMP_MAKE_PHONONS")?strtod(getenv("G4CMP_MAKE_PHONONS"),0):1.),
    genCharges(getenv("G4CMP_MAKE_CHARGES")?strtod(getenv("G4CMP_MAKE_CHARGES"),0):1.),
    lukeSample(getenv("G4CMP_LUKE_SAMPLE")?strtod(getenv("G4CMP_LUKE_SAMPLE"),0):1.),
    downSample(getenv("G4CMP_DOWN_SAMPLE")?strtod(getenv("G4CMP_DOWN_SAMPLE"),0):1.),
    EminPhonons(getenv("G4CMP_EMIN_PHONONS")?strtod(getenv("G4CMP_EMIN_PHONONS"),0)*eV:0.),
    EminCharges(getenv("G4CMP_EMIN_CHARGES")?strtod(getenv("G4CMP_EMIN_CHARGES"),0)*eV:0.),
    useKVsolver(getenv("G4CMP_USE_KVSOLVER")?atoi(getenv("G4CMP_USE_KVSOLVER")):0),
    fanoEnabled(getenv("G4CMP_FANO_ENABLED")?atoi(getenv("G4CMP_FANO_ENABLED")):1),
    chargeCloud(getenv("G4CMP_CHARGE_CLOUD")?atoi(getenv("G4CMP_CHARGE_CLOUD")):0),
    messenger(new G4CMPConfigMessenger(this)) {
  fPhysicsModelID = G4PhysicsModelCatalog::Register("G4CMP process");
}

G4CMPConfigManager::~G4CMPConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void G4CMPConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}
