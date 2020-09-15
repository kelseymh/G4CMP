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
// 20190404  Change "point" to "point3d" to make way for 2D interpolator.
// 20190508  Move some 2D/3D common features to new base class
// 20190630  Have MatInv() return error (false), catch up calling chain.
// 20190923  Add constructor with neighbors table, use with Clone().
// 20200907  Add BuildTInverse() function to precompute invT for Cart2Bary().
//		Add "quiet" argument to MatInv to suppress warnings.
// 20200908  Replace four-arg ctor and UseMesh() with copy constructor.
// 20200914  Include gradient precalculation in BuildTInverse action.

#ifndef G4CMPTriLinearInterp_h 
#define G4CMPTriLinearInterp_h 

#include "G4CMPVMeshInterpolator.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <array>

// Convenient abbreviations, available to subclasses and client code
using mat3x3 = std::array<std::array<G4double,3>,3>;
using mat4x3 = std::array<std::array<G4double,3>,4>;

class G4CMPTriLinearInterp : public G4CMPVMeshInterpolator {
public:
  // Uninitialized version; user MUST call UseMesh()
  G4CMPTriLinearInterp() : G4CMPVMeshInterpolator("TRI") {;}

  // Mesh coordinates and values only; uses QHull to generate triangulation
  G4CMPTriLinearInterp(const std::vector<point3d>& xyz,
		       const std::vector<G4double>& v);

  // Mesh points and pre-defined triangulation
  G4CMPTriLinearInterp(const std::vector<point3d>& xyz,
		       const std::vector<G4double>& v,
		       const std::vector<tetra3d>& tetra);

  // Cloning function to allow making type-matched copies
  virtual G4CMPVMeshInterpolator* Clone() const {
    return new G4CMPTriLinearInterp(*this);
  }

  G4CMPTriLinearInterp(const G4CMPTriLinearInterp& rhs);

  // User initialization or re-initialization
  void UseMesh(const std::vector<point3d>& xyz, const std::vector<G4double>& v);

  void UseMesh(const std::vector<point3d>& xyz, const std::vector<G4double>& v,
	       const std::vector<tetra3d>& tetra);

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  G4double GetValue(const G4double pos[], G4bool quiet=false) const;
  G4ThreeVector GetGrad(const G4double pos[], G4bool quiet=false) const;

  void SavePoints(const G4String& fname) const;
  void SaveTetra(const G4String& fname) const;

protected:
  void FillGradients();		// Compute gradient (field) at each tetrahedron

private:
  std::vector<point3d> X;
  std::vector<tetra3d> Tetrahedra;
  std::vector<tetra3d> Neighbors;
  std::vector<mat3x3> TInverse;		// Matrix for barycenter calculation
  std::vector<mat4x3> TExtend;		// Matrix for gradient calculation
  std::vector<G4bool> TInvGood;		// Flags for noninvertible matrix

  mutable std::map<G4int,G4int> qhull2x;	// Used by QHull for meshing

  // Lists of tetrahedra with shared vertices, for generating neighbors table
  std::vector<tetra3d> Tetra012;	// Duplicate tetrahedra lists
  std::vector<tetra3d> Tetra013;	// Sorted on vertex triplets
  std::vector<tetra3d> Tetra023;
  std::vector<tetra3d> Tetra123;

  void BuildTetraMesh();	// Builds mesh from pre-initialized 'X' array
  void FillNeighbors();		// Generate Neighbors table from tetrahedra
  void FillTInverse();		// Compute inverse matrices for Cart2Bary()

  // Function pointer for comparison operator to use search for facets
  using TetraComp = G4bool(*)(const tetra3d&, const tetra3d&);

  G4int FindNeighbor(const std::array<G4int,3>& facet, G4int skip) const;
  G4int FindTetraID(const std::vector<tetra3d>& tetras,
		    const tetra3d& wildTetra, G4int skip,
		    TetraComp tLess) const;
  G4int FirstInteriorTetra();	// Lowest tetra index with all facets shared

  void FindTetrahedron(const G4double point[3], G4double bary[4],
		       G4bool quiet=false) const;
  G4int FindPointID(const std::vector<G4double>& point, const G4int id) const;

  G4bool Cart2Bary(const G4double point[3], G4double bary[4]) const;
  G4bool BuildT4x3(size_t itet, mat4x3& ET) const;

  G4bool MatInv(const mat3x3& matrix, mat3x3& result, G4bool quiet=false) const;
  G4double BaryNorm(G4double bary[4]) const;
  G4double Det3(const mat3x3& matrix) const;

  // Dump tetrahedron information (neighbors and vertices)
  void PrintTetra(std::ostream& os, G4int iTetra) const;
};

#endif	/* G4CMPTriLinearInterp */
