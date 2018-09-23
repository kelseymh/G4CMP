/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170525  Drop unnecessary empty destructor ("rule of five" pattern)
// 20180525  Add "quiet" flag to suppress "outside of hull" messages
// 20180904  Add constructor to directly load mesh definitions

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

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  G4double GetValue(const G4double pos[3], G4bool quiet=false) const;
  G4ThreeVector GetGrad(const G4double pos[3], G4bool quiet=false) const;
  
private:
  std::vector<point > X;
  std::vector<G4double> V;
  std::vector<std::array<G4int, 4> > Tetrahedra;
  std::vector<std::array<G4int, 4> > Neighbors;
  mutable std::map<G4int,G4int> qhull2x;
  mutable G4int TetraIdx;
  mutable G4ThreeVector cachedGrad;
  mutable G4bool staleCache;

  void BuildTetraMesh();	// Builds mesh from pre-initialized 'X' array
  void FillNeighbors();		// Generate Neighbors table from tetrahedra

  G4int FindNeighbor(const std::array<G4int,3>& facet, G4int skip) const;
  G4int FindTetraID(const std::array<G4int,4>& wildTetra, G4int skip) const;

  void FindTetrahedron(const G4double point[4], G4double bary[4],
		       G4bool quiet=false) const;
  G4int FindPointID(const std::vector<G4double>& point, const G4int id) const;

  void Cart2Bary(const G4double point[4], G4double bary[4]) const;
  void BuildT4x3(G4double ET[4][3]) const;
  G4double Det3(const G4double matrix[3][3]) const;
  void MatInv(const G4double matrix[3][3], G4double result[3][3]) const;
};

#endif
