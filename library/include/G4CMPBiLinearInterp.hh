/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20190404  Adapted from TriLinearInterp for use with 2D triangular mesh

#ifndef G4CMPBiLinearInterp_h 
#define G4CMPBiLinearInterp_h 

#include "globals.hh"
#include "G4TwoVector.hh"
#include <vector>
#include <map>
#include <array>

using point2d = std::array<G4double,2>;


class G4CMPBiLinearInterp {
public:
  // Uninitialized version; user MUST call UseMesh()
  G4CMPBiLinearInterp() : TetraIdx(0), staleCache(true) {;}

  // Mesh points and pre-defined triangulation
  G4CMPBiLinearInterp(const std::vector<point2d>& xy,
		      const std::vector<G4double>& v,
		      const std::vector<std::array<G4int,3> >& tetra);

  // User initialization or re-initialization
  void UseMesh(const std::vector<point2d>& xy, const std::vector<G4double>& v,
	       const std::vector<std::array<G4int,3> >& tetra);

  // Replace values at mesh points without rebuilding tables
  void UseValues(const std::vector<G4double>& v);

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  G4double GetValue(const G4double pos[2], G4bool quiet=false) const;
  G4TwoVector GetGrad(const G4double pos[2], G4bool quiet=false) const;

  void SavePoints(const G4String& fname) const;
  void SaveTetra(const G4String& fname) const;

private:
  std::vector<point2d > X;
  std::vector<G4double> V;
  std::vector<std::array<G4int,3> > Tetrahedra;	// For 2D, these are triangles!
  std::vector<std::array<G4int,3> > Neighbors;
  mutable G4int TetraIdx;
  mutable G4TwoVector cachedGrad;
  mutable G4bool staleCache;
  G4int TetraStart;				// Start of tetrahedral searches

  std::vector<std::array<G4int,3> > Tetra01;	// Duplicate tetrahedra lists
  std::vector<std::array<G4int,3> > Tetra02;	// Sorted on vertex triplets
  std::vector<std::array<G4int,3> > Tetra12;

  void FillNeighbors();		// Generate Neighbors table from tetrahedra

  // Function pointer for comparison operator to use search for facets
  using TetraComp = G4bool(*)(const std::array<G4int,3>&,
			      const std::array<G4int,3>&);

  G4int FindNeighbor(const std::array<G4int,2>& edge, G4int skip) const;
  G4int FindTetraID(const std::vector<std::array<G4int,3> >& tetras,
		    const std::array<G4int,3>& wildTetra, G4int skip,
		    TetraComp tLess) const;
  G4int FirstInteriorTetra();	// Lowest tetra index with all facets shared

  void FindTetrahedron(const G4double point[2], G4double bary[3],
		       G4bool quiet=false) const;

  void Cart2Bary(const G4double point[2], G4double bary[3]) const;
  G4double BaryNorm(G4double bary[3]) const;
  void BuildT3x2(G4double ET[3][2]) const;
  G4double Det2(const G4double matrix[2][2]) const;
  void MatInv(const G4double matrix[2][2], G4double result[2][2]) const;
};

// SPECIAL:  Provide a way to write out array data directly (not in STL!)

template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& arr) {
  for (const T& ai: arr) os << ai << " ";
  return os;
}

#endif
