/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170525  Drop unnecessary empty destructor ("rule of five" pattern)
// 20180525  Add "quiet" flag to suppress "outside of hull" messages
// 20180904  Add constructor to directly load mesh definitions
// 20180926  Add functions to write points, tetrahedra etc. to files.
//		Add starting index for tetrahedral traversal
// 20190226  Provide accessor to replace potentials at mesh points

#ifndef G4CMPTriLinearInterp_h 
#define G4CMPTriLinearInterp_h 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <array>

using point = std::array<G4double, 3>;

class G4CMPTriLinearInterp {
public:
  // Uninitialized version; user MUST call UseMesh()
  G4CMPTriLinearInterp() : TetraIdx(0), staleCache(true) {;}

  // Mesh coordinates and values only; uses QHull to generate triangulation
  G4CMPTriLinearInterp(const std::vector<point>& xyz,
		       const std::vector<G4double>& v);

  // Mesh points and pre-defined triangulation
  G4CMPTriLinearInterp(const std::vector<point>& xyz,
		       const std::vector<G4double>& v,
		       const std::vector<std::array<G4int,4> >& tetra);

  // User initialization or re-initialization
  void UseMesh(const std::vector<point>& xyz, const std::vector<G4double>& v);

  void UseMesh(const std::vector<point>& xyz, const std::vector<G4double>& v,
	       const std::vector<std::array<G4int,4> >& tetra);

  // Replace values at mesh points without rebuilding tables
  void UseValues(const std::vector<G4double>& v);

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  G4double GetValue(const G4double pos[3], G4bool quiet=false) const;
  G4ThreeVector GetGrad(const G4double pos[3], G4bool quiet=false) const;

  void SavePoints(const G4String& fname) const;
  void SaveTetra(const G4String& fname) const;

private:
  std::vector<point > X;
  std::vector<G4double> V;
  std::vector<std::array<G4int,4> > Tetrahedra;
  std::vector<std::array<G4int,4> > Neighbors;
  mutable std::map<G4int,G4int> qhull2x;
  mutable G4int TetraIdx;
  mutable G4ThreeVector cachedGrad;
  mutable G4bool staleCache;
  G4int TetraStart;				// Start of tetrahedral searches

  std::vector<std::array<G4int,4> > Tetra012;	// Duplicate tetrahedra lists
  std::vector<std::array<G4int,4> > Tetra013;	// Sorted on vertex triplets
  std::vector<std::array<G4int,4> > Tetra023;
  std::vector<std::array<G4int,4> > Tetra123;

  void BuildTetraMesh();	// Builds mesh from pre-initialized 'X' array
  void FillNeighbors();		// Generate Neighbors table from tetrahedra

  // Function pointer for comparison operator to use search for facets
  using TetraComp = G4bool(*)(const std::array<G4int,4>&,
			      const std::array<G4int,4>&);

  G4int FindNeighbor(const std::array<G4int,3>& facet, G4int skip) const;
  G4int FindTetraID(const std::vector<std::array<G4int,4> >& tetras,
		    const std::array<G4int,4>& wildTetra, G4int skip,
		    TetraComp tLess) const;
  G4int FirstInteriorTetra();	// Lowest tetra index with all facets shared

  void FindTetrahedron(const G4double point[4], G4double bary[4],
		       G4bool quiet=false) const;
  G4int FindPointID(const std::vector<G4double>& point, const G4int id) const;

  void Cart2Bary(const G4double point[4], G4double bary[4]) const;
  G4double BaryNorm(G4double bary[4]) const;
  void BuildT4x3(G4double ET[4][3]) const;
  G4double Det3(const G4double matrix[3][3]) const;
  void MatInv(const G4double matrix[3][3], G4double result[3][3]) const;
};

// SPECIAL:  Provide a way to write out array data directly (not in STL!)

template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& arr) {
  for (const T& ai: arr) os << ai << " ";
  return os;
}

#endif
