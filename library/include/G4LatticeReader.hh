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

#ifndef G4LatticeReader_h
#define G4LatticeReader_h 1

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>

class G4LatticeLogical;

class G4LatticeReader {
public:
  G4LatticeReader(G4int vb=0);
  ~G4LatticeReader();

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  G4LatticeLogical* MakeLattice(const G4String& filepath);

protected:
  G4bool OpenFile(const G4String& filepath);
  G4bool ProcessToken();
  G4bool ProcessValue(const G4String& name);	// Numerical parameters
  G4bool ProcessConstants();			// Four dynamical constants
  G4bool ProcessBasisVector();			// Crystal basis vector
  G4bool ProcessMassTensor();			// Electron mass tensor
  G4bool ProcessEulerAngles(const G4String& name);	// Drift directions
  G4bool ProcessMap();				// Velocity magnitudes file
  G4bool ProcessNMap();				// Direction vectors file
  G4bool ReadMapInfo();				// Get map file parameters
  G4bool SkipComments();			// Everything after '#'
  void CloseFile();

private:
  G4int verboseLevel;		// For reporting progress, also use G4VERBOSE

  std::ifstream* psLatfile;	// Configuration file being read
  G4LatticeLogical* pLattice;	// Lattice under construction (not owned)

  G4String fMapPath;		// Path to config file to find velocity maps
  G4String fToken;		// Reusable buffers for reading file
  G4double fValue;		// ... floating point data values
  G4String fMap, fsPol;		// ... map filename and polarization code
  G4int    fPol, fNX, fNY;	// ... map binning in each direction
  G4RotationMatrix fMatrix;	// ... 3x3 matrix for mass, drift valleys
  G4ThreeVector f3Vec;		// ... three-vector for mass, crystal bases
  G4int fLastBasis;		// ... index counter for basis vectors

  const G4String fDataDir;	// Directory path ($G4LATTICEDATA)
  const G4double mElectron;	// Electron mass in kilograms
};

#endif	/* G4LatticeReader_h */
