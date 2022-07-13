/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/src/G4LatticeManager.cc
/// \brief Implementation of the G4LatticeManager class
//
// $Id$
//
// 20131113  Delete lattices in (new) registry, not in lookup maps
// 20140319  Dump lattices on load with verbosity; don't double register!
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20140412  Use const volumes and materials for registration
// 20141008  Change to global singleton; must be shared across worker threads
// 20160111  Remove G4 version checking. We hard depend on 10.2+
// 20160615  Set name of G4LatticeLogical to match config directory
// 20160826  Get default verbosity from envvar
// 20170527  Drop unnecessary <fstream>
// 20170817  Increase verbosity cut on informational messages
// 20170928  Replace "polarizationState" with "mode"

#include "G4LatticeManager.hh"
#include "G4CMPConfigManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4LatticeReader.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

G4LatticeManager* G4LatticeManager::fLM = 0;

#include "G4AutoLock.hh"
namespace {
  G4Mutex latMutex = G4MUTEX_INITIALIZER;     // For thread protection
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeManager::G4LatticeManager()
  : verboseLevel(G4CMPConfigManager::GetVerboseLevel()) {
  Clear();
}

G4LatticeManager::~G4LatticeManager() {
  Reset();			// Deletes all lattices
}

// Delete all registered lattices and clear entries from lookup tables

void G4LatticeManager::Reset() {
  for (LatticeLogReg::iterator lm=fLLattices.begin();
       lm != fLLattices.end(); ++lm) {
    delete (*lm);
  }

  for (LatticePhyReg::iterator pm=fPLattices.begin();
       pm != fPLattices.end(); ++pm) {
    delete (*pm);
  }

  Clear();
}

// Remove entries without deletion (for begin-job and end-job initializing)

void G4LatticeManager::Clear() {
  fPLatticeList.clear();
  fPLattices.clear();

  fLLatticeList.clear();
  fLLattices.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeManager* G4LatticeManager::GetLatticeManager() {
  G4AutoLock latLock(&latMutex);      // Protect before changing pointer

  // if no lattice manager exists, create one.
  if (!fLM) fLM = new G4LatticeManager();
  return fLM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Associate logical lattice with material

G4bool G4LatticeManager::RegisterLattice(const G4Material* Mat,
					 G4LatticeLogical* Lat) {
  if (!Mat || !Lat) return false;	// Don't register null pointers

  if (fLLatticeList.find(Mat) != fLLatticeList.end()) {
    G4cerr << "WARNING:  Replacing existing lattice for " << Mat->GetName()
	   << G4endl;
  }

  G4AutoLock latLock(&latMutex);      // Protect before changing registry

  fLLattices.insert(Lat);		// Take ownership in registry
  fLLatticeList[Mat] = Lat;

  if (verboseLevel) {
    G4cout << "G4LatticeManager::RegisterLattice: "
	   << " Total number of logical lattices: " << fLLatticeList.size()
	   << " (" << fLLattices.size() << " unique)" << G4endl;
  }

  return true; 
}

// Construct logical lattice for material from config file

G4LatticeLogical* G4LatticeManager::LoadLattice(const G4Material* Mat,
						const G4String& latDir) {
  if (verboseLevel) {
    G4cout << "G4LatticeManager::LoadLattice material " << Mat->GetName()
	   << " " << latDir << G4endl;
  }
		      
  G4LatticeReader latReader(verboseLevel);
  G4LatticeLogical* newLat = latReader.MakeLattice(latDir+"/config.txt");
  if (newLat) {
    newLat->SetDensity(Mat->GetDensity());
    newLat->Initialize(latDir);
    if (verboseLevel>1)
      G4cout << " Created newLat " << newLat << "\n" << *newLat << G4endl;

    RegisterLattice(Mat, newLat);
  } else {
    G4cerr << "ERROR creating " << latDir << " lattice for material "
	   << Mat->GetName() << G4endl;
  }

  return newLat;
}

// Combine loading and registration (Material extracted from volume)

G4LatticePhysical* G4LatticeManager::LoadLattice(const G4VPhysicalVolume* Vol,
						const G4String& latDir) {
  if (verboseLevel) {
    G4cout << "G4LatticeManager::LoadLattice volume " << Vol->GetName()
	   << " " << latDir << G4endl;
  }
		      
  const G4Material* theMat = Vol->GetLogicalVolume()->GetMaterial();

  // Create and register the logical lattice, then the physical lattice
  G4LatticeLogical* lLattice = LoadLattice(theMat, latDir);
  if (!lLattice) return 0;

  G4LatticePhysical* pLattice = new G4LatticePhysical(lLattice);
  if (pLattice) RegisterLattice(Vol, pLattice);

  if (verboseLevel>1)
    G4cout << " Created pLattice " << pLattice << "\n" << *pLattice << G4endl;

  return pLattice;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Associate physical (oriented) lattice with physical volume

G4bool G4LatticeManager::RegisterLattice(const G4VPhysicalVolume* Vol,
					 G4LatticePhysical* Lat) {
  if (!Vol || !Lat) return false;	// Don't register null pointers

  G4AutoLock latLock(&latMutex);      // Protect before changing pointer

  // SPECIAL: Register first lattice with a null volume to act as default
  if (fPLatticeList.empty()) fPLatticeList[0] = Lat;

  if (fPLatticeList.find(Vol) != fPLatticeList.end()) {
    G4cerr << "WARNING:  Replacing existing lattice for " << Vol->GetName()
	   << G4endl;
  }

  fPLattices.insert(Lat);
  fPLatticeList[Vol] = Lat;

  if (verboseLevel) {
    G4cout << "G4LatticeManager::RegisterLattice: "
	   << " Total number of physical lattices: " << fPLatticeList.size()-1
	   << " (" << fPLattices.size() << " unique)" << G4endl;
  }

  return true; 
}

G4bool G4LatticeManager::RegisterLattice(const G4VPhysicalVolume* Vol,
					 G4LatticeLogical* LLat) {
  if (!Vol || !LLat) return false;	// Don't register null pointers

  // Make sure logical lattice is registered for material
  RegisterLattice(Vol->GetLogicalVolume()->GetMaterial(), LLat);

  // Create and register new physical lattice to go with volume
  return RegisterLattice(Vol, new G4LatticePhysical(LLat));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Returns a pointer to the LatticeLogical associated with material

G4LatticeLogical* G4LatticeManager::GetLattice(const G4Material* Mat) const {
  LatticeMatMap::const_iterator latFind = fLLatticeList.find(Mat);
  if (latFind != fLLatticeList.end()) {
    if (verboseLevel>2)
      G4cout << "G4LatticeManager::GetLattice found " << latFind->second
	     << " for " << (Mat?Mat->GetName():"NULL") << "." << G4endl;
    return latFind->second;
  }


  if (verboseLevel) 
    G4cerr << "G4LatticeManager:: Found no matching lattices for "
	   << (Mat?Mat->GetName():"NULL") << "." << G4endl;

  return 0;			// No lattice associated with volume
}

// Returns a pointer to the LatticePhysical associated with volume
// NOTE:  Passing Vol==0 will return the default lattice

G4LatticePhysical* 
G4LatticeManager::GetLattice(const G4VPhysicalVolume* Vol) const {
  LatticeVolMap::const_iterator latFind = fPLatticeList.find(Vol);
  if (latFind != fPLatticeList.end()) {
    if (verboseLevel>2)
      G4cout << "G4LatticeManager::GetLattice found " << latFind->second
	     << " for " << (Vol?Vol->GetName():"default") << "." << G4endl;
    return latFind->second;
  }

  if (verboseLevel) 
    G4cerr << "G4LatticeManager::GetLattice found no matching lattices for "
	   << (Vol?Vol->GetName():"default") << "." << G4endl;

  return 0;			// No lattice associated with volume
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Return true if volume Vol has a physical lattice

G4bool G4LatticeManager::HasLattice(const G4VPhysicalVolume* Vol) const {
  return (fPLatticeList.find(Vol) != fPLatticeList.end());
}

// Return true if material Mat has a logical lattice

G4bool G4LatticeManager::HasLattice(const G4Material* Mat) const {
  return (fLLatticeList.find(Mat) != fLLatticeList.end());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and mode(0=LON, 1=FT, 2=ST), 
//returns phonon velocity in m/s

G4double G4LatticeManager::MapKtoV(const G4VPhysicalVolume* Vol, G4int mode,
				   const G4ThreeVector & k) const {
  G4LatticePhysical* theLattice = GetLattice(Vol);
  if (verboseLevel>2)
    G4cout << "G4LatticeManager::MapKtoV using lattice " << theLattice
	   << G4endl;

  // If no lattice available, use generic "speed of sound"
  return theLattice ? theLattice->MapKtoV(mode, k) : 300.*m/s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Given the phonon wave vector k, phonon physical volume Vol 
// and mode(0=LON, 1=FT, 2=ST), 
// returns phonon propagation direction as dimensionless unit vector

G4ThreeVector 
G4LatticeManager::MapKtoVDir(const G4VPhysicalVolume* Vol, G4int mode,
			     const G4ThreeVector & k) const {
  G4LatticePhysical* theLattice = GetLattice(Vol);
  if (verboseLevel>2)
    G4cout << "G4LatticeManager::MapKtoVDir using lattice " << theLattice
	   << G4endl;

  // If no lattice available, propagate along input wavevector
  return theLattice ? theLattice->MapKtoVDir(mode, k) : k.unit();
}
