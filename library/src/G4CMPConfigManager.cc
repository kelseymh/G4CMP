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
// 20180801  G4CMP-143:  Change IV rate from bool to str, Edelweiss->Quadratic
// 20190711  G4CMP-158:  Add functions to select NIEL yield functions
// 20191014  G4CMP-179:  Drop sampling of anharmonic decay (downconversion)
// 20200211  G4CMP-191:  Add version identification from .g4cmp-version
// 20200331  G4CMP-195:  Add charge trapping mean free paths
// 20200331  G4CMP-196:  Add impact ionization mean free path
// 20200426  G4CMP-196: Change "impact ionization" to "trap ionization"
// 20200501  G4CMP-196: Change trap-ionization MFP names, "eTrap" -> "DTrap",
//		"hTrap" -> "ATrap".
// 20200504  G4CMP-195:  Reduce length of charge-trapping parameter names
// 20200530  G4CMP-202:  Provide separate master and worker instances
// 20200614  G4CMP-211:  Add functionality to print settings
// 20200614  G4CMP-210:  Add missing initializers to copy constructor
// 20210303  G4CMP-243:  Add parameter to set step length for merging hits
// 20210910  G4CMP-272:  Add parameter to set number of downsampled Luke phonons
// 20220921  G4CMP-319:  Add temperature setting for use with QP sensors.
// 20221014  G4CMP-334:  Add maxLukePhonons to printout; show macro commands

#include "G4CMPConfigManager.hh"
#include "G4CMPConfigMessenger.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4CMPLindhardNIEL.hh"
#include "G4VNIELPartition.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <float.h>
#include <stdlib.h>


// Singleton Initializers for master and worker threads

G4CMPConfigManager* G4CMPConfigManager::Instance() {
  static const G4CMPConfigManager* masterInstance = 0;
  static G4ThreadLocal G4CMPConfigManager* theInstance = 0;

  if (!theInstance) {
    if (!G4Threading::IsWorkerThread()) {	// Master or sequential
      theInstance = new G4CMPConfigManager;
      masterInstance = theInstance;
    } else {					// Workers copy from master
      theInstance = new G4CMPConfigManager(*masterInstance);
    }
  }

  return theInstance;
}

// Object constructor

G4CMPConfigManager::G4CMPConfigManager()
  : verbose(getenv("G4CMP_DEBUG")?atoi(getenv("G4CMP_DEBUG")):0),
    ehBounces(getenv("G4CMP_EH_BOUNCES")?atoi(getenv("G4CMP_EH_BOUNCES")):1),
    pBounces(getenv("G4CMP_PHON_BOUNCES")?atoi(getenv("G4CMP_PHON_BOUNCES")):100),
    maxLukePhonons(getenv("G4MP_MAX_LUKE")?atoi(getenv("G4MP_MAX_LUKE")):-1),
    LatticeDir(getenv("G4LATTICEDATA")?getenv("G4LATTICEDATA"):"./CrystalMaps"),
    IVRateModel(getenv("G4CMP_IV_RATE_MODEL")?getenv("G4CMP_IV_RATE_MODEL"):"Quadratic"),
    eTrapMFP(getenv("G4CMP_ETRAPPING_MFP")?strtod(getenv("G4CMP_ETRAPPING_MFP"),0)*mm:DBL_MAX),
    hTrapMFP(getenv("G4CMP_HTRAPPING_MFP")?strtod(getenv("G4CMP_HTRAPPING_MFP"),0)*mm:DBL_MAX),
    eDTrapIonMFP(getenv("G4CMP_EDTRAPION_MFP")?strtod(getenv("G4CMP_EDTRAPION_MFP"),0)*mm:DBL_MAX),
    eATrapIonMFP(getenv("G4CMP_EATRAPION_MFP")?strtod(getenv("G4CMP_EATRAPION_MFP"),0)*mm:DBL_MAX),
    hDTrapIonMFP(getenv("G4CMP_HDTRAPION_MFP")?strtod(getenv("G4CMP_HDTRAPION_MFP"),0)*mm:DBL_MAX),
    hATrapIonMFP(getenv("G4CMP_HATRAPION_MFP")?strtod(getenv("G4CMP_HATRAPION_MFP"),0)*mm:DBL_MAX),
    temperature(getenv("G4CMP_TEMPERATURE")?strtod(getenv("G4CMP_TEMPERATURE"),0)*kelvin:0.),
    clearance(getenv("G4CMP_CLEARANCE")?strtod(getenv("G4CMP_CLEARANCE"),0)*mm:1e-6*mm),
    stepScale(getenv("G4CMP_MIN_STEP")?strtod(getenv("G4CMP_MIN_STEP"),0):-1.),
    sampleEnergy(getenv("G4CMP_SAMPLE_ENERGY")?strtod(getenv("G4CMP_SAMPLE_ENERGY"),0):-1.),
    genPhonons(getenv("G4CMP_MAKE_PHONONS")?strtod(getenv("G4CMP_MAKE_PHONONS"),0):1.),
    genCharges(getenv("G4CMP_MAKE_CHARGES")?strtod(getenv("G4CMP_MAKE_CHARGES"),0):1.),
    lukeSample(getenv("G4CMP_LUKE_SAMPLE")?strtod(getenv("G4CMP_LUKE_SAMPLE"),0):1.),
    combineSteps(getenv("G4CMP_COMBINE_STEPLEN")?strtod(getenv("G4CMP_COMBINE_STEPLEN"),0):0.),
    EminPhonons(getenv("G4CMP_EMIN_PHONONS")?strtod(getenv("G4CMP_EMIN_PHONONS"),0)*eV:0.),
    EminCharges(getenv("G4CMP_EMIN_CHARGES")?strtod(getenv("G4CMP_EMIN_CHARGES"),0)*eV:0.),
    useKVsolver(getenv("G4CMP_USE_KVSOLVER")?atoi(getenv("G4CMP_USE_KVSOLVER")):0),
    fanoEnabled(getenv("G4CMP_FANO_ENABLED")?atoi(getenv("G4CMP_FANO_ENABLED")):1),
    chargeCloud(getenv("G4CMP_CHARGE_CLOUD")?atoi(getenv("G4CMP_CHARGE_CLOUD")):0),
    nielPartition(0), messenger(new G4CMPConfigMessenger(this)) {
  fPhysicsModelID = G4PhysicsModelCatalog::Register("G4CMP process");

  setVersion();

  if (getenv("G4CMP_NIEL_FUNCTION")) 
    setNIEL(getenv("G4CMP_NIEL_FUNCTION"));
  else 
    setNIEL(new G4CMPLewinSmithNIEL);
}

G4CMPConfigManager::~G4CMPConfigManager() {
  delete messenger; messenger=0;
}

// Duplicate existing (master) instances; don't need to check envvars

G4CMPConfigManager::G4CMPConfigManager(const G4CMPConfigManager& master)
  : verbose(master.verbose), fPhysicsModelID(master.fPhysicsModelID), 
    ehBounces(master.ehBounces), pBounces(master.pBounces),
    maxLukePhonons(master.maxLukePhonons),
    version(master.version), LatticeDir(master.LatticeDir), 
    IVRateModel(master.IVRateModel), eTrapMFP(master.eTrapMFP),
    hTrapMFP(master.hTrapMFP), eDTrapIonMFP(master.eDTrapIonMFP),
    eATrapIonMFP(master.eATrapIonMFP), hDTrapIonMFP(master.hDTrapIonMFP),
    hATrapIonMFP(master.hATrapIonMFP),
    temperature(master.temperature), clearance(master.clearance), 
    stepScale(master.stepScale), sampleEnergy(master.sampleEnergy), 
    genPhonons(master.genPhonons), genCharges(master.genCharges), 
    lukeSample(master.lukeSample), combineSteps(master.combineSteps),
    EminPhonons(master.EminPhonons), EminCharges(master.EminCharges),
    useKVsolver(master.useKVsolver), fanoEnabled(master.fanoEnabled),
    chargeCloud(master.chargeCloud), nielPartition(master.nielPartition),
    messenger(new G4CMPConfigMessenger(this)) {;}


// Trigger rebuild of geometry if parameters change

void G4CMPConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}


// Read version tag at build time from generated .g4cmp-version file

void G4CMPConfigManager::setVersion() {
  G4String dir = getenv("G4CMPINSTALL") ? getenv("G4CMPINSTALL") : ".";

  std::ifstream ver(dir+"/.g4cmp-version");
  if (ver.good()) ver >> version;
  else version = "";
}


// Convert input name string to NIEL partitioning function

void G4CMPConfigManager::setNIEL(G4String name) {
  name.toLower();
  if (name(0,3) == "lin") setNIEL(new G4CMPLindhardNIEL);
  if (name(0,3) == "lew") setNIEL(new G4CMPLewinSmithNIEL);
}

void G4CMPConfigManager::setNIEL(G4VNIELPartition* niel) {
  delete nielPartition;
  nielPartition = niel;
}


// Report configuration setting for diagnostics

void G4CMPConfigManager::printConfig(std::ostream& os) const {
  os << "G4CMPConfigManager for G4CMP Version " << version
     << "\nfPhysicsModelID " << fPhysicsModelID
     << "\n/g4cmp/LatticeData " << LatticeDir << "\t# G4LATTICEDATA"
     << "\n/g4cmp/verbose " << verbose << "\t\t\t\t# G4CMP_DEBUG"
     << "\n/g4cmp/chargeBounces " << ehBounces << "\t\t\t\t# G4CMP_EH_BOUNCES"
     << "\n/g4cmp/phononBounces " << pBounces << "\t\t\t# G4CMP_PHON_BOUNCES"
     << "\n/g4cmp/IVRateModel " << IVRateModel << "\t\t\t# G4CMP_IV_RATE_MODEL"
     << "\n/g4cmp/eTrappingMFP " << eTrapMFP/mm << " mm\t\t# G4CMP_ETRAPPING_MFP"
     << "\n/g4cmp/hTrappingMFP " << hTrapMFP/mm << " mm\t\t# G4CMP_HTRAPPING_MFP"
     << "\n/g4cmp/eDTrapIonizationMFP " << eDTrapIonMFP/mm << " mm\t# G4CMP_EDTRAPION_MFP"
     << "\n/g4cmp/eATrapIonizationMFP " << eATrapIonMFP/mm << " mm\t# G4CMP_EATRAPION_MFP"
     << "\n/g4cmp/hDTrapIonizationMFP " << hDTrapIonMFP/mm << " mm\t# G4CMP_HDTRAPION_MFP"
     << "\n/g4cmp/hATrapIonizationMFP " << hATrapIonMFP/mm << " mm\t# G4CMP_HATRAPION_MFP"
     << "\n/g4cmp/temperature " << temperature/kelvin << " K\t\t\t\t# G4CMP_TEMPERATURE"
     << "\n/g4cmp/clearance " << clearance/mm << " mm\t\t\t# G4CMP_CLEARANCE"
     << "\n/g4cmp/minimumStep " << stepScale << "\t\t\t\t# G4CMP_MIN_STEP"
     << "\n/g4cmp/samplingEnergy " << sampleEnergy << "\t\t\t# G4CMP_SAMPLE_ENERGY"
     << "\n/g4cmp/producePhonons " << genPhonons << "\t\t\t\t# G4CMP_MAKE_PHONONS"
     << "\n/g4cmp/produceCharges " << genCharges << "\t\t\t\t# G4CMP_MAKE_CHARGES"
     << "\n/g4cmp/sampleLuke " << lukeSample << "\t\t\t\t# G4CMP_LUKE_SAMPLE"
     << "\n/g4cmp/maxLukePhonons " << maxLukePhonons << "\t\t\t# G4CMP_MAX_LUKE"
     << "\n/g4cmp/combiningStepLength " << combineSteps/mm << " mm\t\t\t# G4CMP_COMBINE_STEPLEN"
     << "\n/g4cmp/minEPhonons " << EminPhonons/eV << " eV\t\t\t\t# G4CMP_EMIN_PHONONS"
     << "\n/g4cmp/minECharges " << EminCharges/eV << " eV\t\t\t\t# G4CMP_EMIN_CHARGES"
     << "\n/g4cmp/useKVsolver " << useKVsolver << "\t\t\t\t# G4CMP_USE_KVSOLVER"
     << "\n/g4cmp/enableFanoStatistics " << fanoEnabled << "\t\t\t# G4CMP_FANO_ENABLED"
     << "\n/g4cmp/createChargeCloud " << chargeCloud << "\t\t\t# G4CMP_CHARGE_CLOUD"
     << "\n/g4cmp/NIELPartition "
     << (nielPartition ? typeid(*nielPartition).name() : "---")
     << "\t# G4CMP_NIEL_FUNCTION "
     << std::endl;
}
