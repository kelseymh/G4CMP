/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// G4CMPVMeshInterpolator:  Base interface class for BiLinear and TriLinear
// interpolators, used by G4CMPMeshElectricField.  This base class allows
// use of a single pointer with a type selected at runtime for the field
// model.  The actual interpolation functionality is implemented entirely
// in the concrete subclasses.
//
// 20200908  Add operator<<() to print matrices (array of array)

#ifndef G4CMPVMeshInterpolator_h 
#define G4CMPVMeshInterpolator_h 

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include <array>
#include <vector>

// Convenient abbreviations, available to subclasses and client code
using point2d = std::array<G4double,2>;
using point3d = std::array<G4double,3>;
using tetra2d = std::array<G4int,3>;
using tetra3d = std::array<G4int,4>;


class G4CMPVMeshInterpolator {
protected:
  // This class CANNOT be instantiated directly!
  G4CMPVMeshInterpolator(const G4String& prefix)
    : TetraIdx(-1), staleCache(true), TetraStart(-1), savePrefix(prefix) {;}

public:
  virtual ~G4CMPVMeshInterpolator() {;}

  // Subclasses MUST implement this function to return duplicate of self
  virtual G4CMPVMeshInterpolator* Clone() const = 0;

public:
  // Replace values at mesh points without rebuilding tables
  void UseValues(const std::vector<G4double>& v);

  // Subclasses MUST implement these functions for their dimensionality

  // Replace existing mesh vectors and tetrahedra table
  // NOTE: Both 2D and 3D versions are given, subclasses should implement one
  void UseMesh(const std::vector<point3d>& /*xyz*/,
	       const std::vector<G4double>& /*v*/,
	       const std::vector<tetra3d>& /*tetra*/) {;}

  void UseMesh(const std::vector<point2d>& /*xy*/,
	       const std::vector<G4double>& /*v*/,
	       const std::vector<tetra2d>& /*tetra*/) {;}

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  virtual G4double GetValue(const G4double pos[], G4bool quiet=false) const = 0;
  virtual G4ThreeVector GetGrad(const G4double pos[], G4bool quiet=false) const = 0;

  // Write out mesh coordinates and tetrahedra table to text files
  virtual void SavePoints(const G4String& fname) const = 0;
  virtual void SaveTetra(const G4String& fname) const = 0;

protected:		// Data members available to subclasses directly
  std::vector<G4double> V;		// Values at mesh points
  // NOTE: Subclasses must define dimensional mesh coords and tetrahera

  mutable G4int TetraIdx;		// Last tetrahedral index used
  mutable G4bool staleCache;		// Flag if cache must be discarded
  mutable G4ThreeVector cachedGrad;
  G4int TetraStart;			// Start of tetrahedral searches

  G4String savePrefix;			// for use in debugging, SaveXxx()
};

// SPECIAL:  Provide a way to write out array/matrix data directly (not in STL!)

template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& arr) {
  for (const T& ai: arr) os << ai << " ";
  return os;
}

template <typename T, size_t M, size_t N>
inline std::ostream& 
operator<<(std::ostream& os, const std::array<std::array<T,N>,M>& mat) {
  for (const auto& ai: mat) os << " " << ai << "\n";
  return os;
}

#endif	/* G4CMPVMeshInterpolator_h */
