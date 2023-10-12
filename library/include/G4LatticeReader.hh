/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4LatticeReader.hh
/// \brief Definition of the G4LatticeReader class
//
// NOTE:  This reader class for logical lattices should be moved to
//	  materials/ after the 10.0 release (and this comment removed).
// $Id$
//
// 20131115  Move ctor, dtor implementations to .cc file.
// 20140218  Add support for charge-carrier functionality
// 20151211  Change fDataDir from static to not static.
// 20160517  Add support to set crystal basis vectors.
// 20160615  Add support to set elasticity tensor.
// 20160630  Drop loading of K-Vg lookup table files
// 20160701  Withdraw seting basis vectors, set crystal symmetry instead
// 20160727  Add functions to handle processing units from config file,
//		define additional units for solid state physics use
// 20170525  Implement 'rule of five' with default copy/move semantics
// 20170810  Add utility function to process list of values with unit.
// 20190704  Add utility function to process string/name argument

#ifndef G4LatticeReader_h
#define G4LatticeReader_h 1

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>
#include <vector>

class G4LatticeLogical;

class G4LatticeReader {
public:
  G4LatticeReader(G4int vb=0);
  ~G4LatticeReader();

  // Use default copy/move semantics
  G4LatticeReader(const G4LatticeReader&) = default;
  G4LatticeReader(G4LatticeReader&&) = default;
  G4LatticeReader& operator=(const G4LatticeReader&) = default;
  G4LatticeReader& operator=(G4LatticeReader&&) = default;

  // Configuration actions
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  G4LatticeLogical* MakeLattice(const G4String& filepath);

protected:
  void DefineUnits();		// Create time^3 and time^4 units for rates

  G4bool OpenFile(const G4String& filepath);
  void CloseFile();

  G4bool ProcessToken();
  G4bool ProcessValue(const G4String& name);	// Single numerical parameter
  G4bool ProcessList(const G4String& unitcat);	// List of parameters with unit
  G4bool ProcessString(const G4String& name);	// Single string parameter

  G4bool ProcessConstants();			// Four dynamical constants
  G4bool ProcessMassTensor();			// Electron mass tensor
  G4bool ProcessCrystalGroup(const G4String& name);	// Symmetry, spacing
  G4bool ProcessDebyeLevel();			// Frequency or temperature
  G4bool ProcessStiffness();			// Elasticity matrix element
  G4bool ProcessEulerAngles(const G4String& name);	// Drift directions
  G4bool ProcessDeformation();			// IV deformation potentials
  G4bool ProcessThresholds();			// IV energy thresholds
  G4bool SkipComments();			// Everything after '#'

  // Read expected dimensions for value from file, return scale factor
  // NOTE: String from file may have leading "/" for inverse units
  // Input argument "unitcat" may be comma-delimited list of categories
  G4double ProcessUnits(const G4String& unitcat);
  G4double ProcessUnits(const G4String& unit, const G4String& unitcat);

private:
  G4int verboseLevel;		// For reporting progress, also use G4VERBOSE

  std::ifstream* psLatfile;	// Configuration file being read
  G4LatticeLogical* pLattice;	// Lattice under construction (not owned)

  G4String fToken;		// Reusable buffers for reading file
  G4double fValue;		// ... floating point data value
  std::vector<G4double> fList;	// ... list of floating point values
  G4RotationMatrix fMatrix;	// ... 3x3 matrix for mass, drift valleys
  G4ThreeVector f3Vec;		// ... three-vector for mass
  G4double fUnits;		// ... dimensional unit scale factor
  G4String fUnitName;		// ... unit string from reading file
  G4String fUnitCat;		// ... G4UnitsCategory of dimensions

  const G4String fDataDir;	// Directory path ($G4LATTICEDATA)
  const G4double mElectron;	// Electron mass in kilograms
};

#endif	/* G4LatticeReader_h */
