/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/include/G4LatticeManager.hh
/// \brief Definition of the G4LatticeManager class
//
// $Id$
//
// 20131113  Add registry to carry unique lattice pointers, for EOJ deletion
// 20131115  Drop lattice counters, not used anywhere
// 20140412  Use const volumes and materials for registration
// 20141008  Change to global singleton; must be shared across worker threads

#ifndef G4LatticeManager_h
#define G4LatticeManager_h 1


#include "G4ThreeVector.hh"
#include <map>
#include <set>

class G4LatticeLogical;
class G4LatticePhysical;
class G4Material;
class G4VPhysicalVolume;


class G4LatticeManager {
private:
  static G4LatticeManager* fLM;		// Singleton

public:
  static G4LatticeManager* GetLatticeManager(); 
  static G4LatticeManager* Instance() { return GetLatticeManager(); }

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  void Reset();		// Remove and delete all registered lattices

  // Users may register physical or logical lattices with volumes
  G4bool RegisterLattice(const G4VPhysicalVolume*, G4LatticePhysical*);
  G4bool RegisterLattice(const G4VPhysicalVolume*, G4LatticeLogical*);

  // Logical lattices are associated with materials
  G4bool RegisterLattice(const G4Material*, G4LatticeLogical*);

  // Logical lattices may be read from <latDir>/config.txt data file
  G4LatticeLogical* LoadLattice(const G4Material*, const G4String& latDir);
  G4LatticeLogical* GetLattice(const G4Material*) const;
  G4bool HasLattice(const G4Material*) const;

  // Combine loading and registration (Material extracted from volume)
  G4LatticePhysical* LoadLattice(const G4VPhysicalVolume*,
				 const G4String& latDir);

  // NOTE:  Passing Vol==0 will return the default lattice
  G4LatticePhysical* GetLattice(const G4VPhysicalVolume*) const;
  G4bool HasLattice(const G4VPhysicalVolume*) const;

  G4double MapKtoV(const G4VPhysicalVolume*, G4int,
		   const G4ThreeVector &) const;

  G4ThreeVector MapKtoVDir(const G4VPhysicalVolume*, G4int,
			   const G4ThreeVector&) const;

protected:
  void Clear();		// Remove entries from lookup tables w/o deletion

protected:
  G4int verboseLevel;		// Allow users to enable diagnostic messages

  typedef std::map<const G4Material*, G4LatticeLogical*> LatticeMatMap;
  typedef std::set<G4LatticeLogical*> LatticeLogReg;

  LatticeLogReg fLLattices;	// Registry of unique lattice pointers
  LatticeMatMap fLLatticeList;

  typedef std::map<const G4VPhysicalVolume*, G4LatticePhysical*> LatticeVolMap;
  typedef std::set<G4LatticePhysical*> LatticePhyReg;

  LatticePhyReg fPLattices;	// Registry of unique lattice pointers
  LatticeVolMap fPLatticeList; 

private:
  G4LatticeManager();
  virtual ~G4LatticeManager();
};

#endif
