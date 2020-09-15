/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20190404  Adapted from TriLinearInterp for use with 2D triangular mesh
// 20190508  Move some 2D/3D common features to new base class
// 20190630  Have MatInv() return error (false), catch up calling chain.
// 20190923  Add constructor with neighbors table, use with Clone().
// 20200907  Add BuildTInverse() function to precompute invT for Cart2Bary().
//             Add "quiet" argument to MatInv to suppress warnings.
// 20200908  Replace four-arg ctor and UseMesh() with copy constructor.
// 20200914  Include gradient precalculation in BuildTInverse action.

#ifndef G4CMPBiLinearInterp_h 
#define G4CMPBiLinearInterp_h 

#include "G4CMPVMeshInterpolator.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <array>

// Convenient abbreviations, available to subclasses and client code
using mat2x2 = std::array<std::array<G4double,2>,2>;
using mat3x2 = std::array<std::array<G4double,2>,3>;
 

class G4CMPBiLinearInterp : public G4CMPVMeshInterpolator {
public:
  // Uninitialized version; user MUST call UseMesh()
  G4CMPBiLinearInterp() : G4CMPVMeshInterpolator("BLI") {;}

  // Mesh points and pre-defined triangulation
  G4CMPBiLinearInterp(const std::vector<point2d>& xy,
		      const std::vector<G4double>& v,
		      const std::vector<tetra2d>& tetra);

  // Allow use of 3D inputs, which get collapsed to 2D internals
  G4CMPBiLinearInterp(const std::vector<point3d>& xyz,
		      const std::vector<G4double>& v,
		      const std::vector<tetra3d>& tetra);

  // Cloning function to allow making type-matched copies
  virtual G4CMPVMeshInterpolator* Clone() const {
    return new G4CMPBiLinearInterp(*this);
  }

  G4CMPBiLinearInterp(const G4CMPBiLinearInterp& rhs);

  // User initialization or re-initialization
  void UseMesh(const std::vector<point2d>& xy, const std::vector<G4double>& v,
	       const std::vector<tetra2d>& tetra);

  void UseMesh(const std::vector<point3d>& xyz, const std::vector<G4double>& v,
	       const std::vector<tetra3d>& tetra);

  // Evaluate mesh at arbitrary location, optionally suppressing errors
  G4double GetValue(const G4double pos[], G4bool quiet=false) const;
  G4ThreeVector GetGrad(const G4double pos[], G4bool quiet=false) const;

  void SavePoints(const G4String& fname) const;
  void SaveTetra(const G4String& fname) const;

private:
  std::vector<point2d> X;
  std::vector<tetra2d> Tetrahedra;	// For 2D, these are triangles!
  std::vector<tetra2d> Neighbors;
  std::vector<mat2x2> TInverse;		// Matrix for barycenter calculation
  std::vector<mat3x2> TExtend;		// Matrix for gradient calculation
  std::vector<G4bool> TInvGood;		// Flags for noninvertible matrix

  std::vector<tetra2d> Tetra01;		// Duplicate tetrahedra lists
  std::vector<tetra2d> Tetra02;		// Sorted on vertex triplets
  std::vector<tetra2d> Tetra12;

  void FillNeighbors();		// Generate Neighbors table from tetrahedra
  void BuildTInverse();		// Compute inverse matrices for Cart2Bary()

  void Compress3DPoints(const std::vector<point3d>& xyz);
  void Compress3DTetras(const std::vector<tetra3d>& tetra);

  // Function pointer for comparison operator to use search for facets
  using TetraComp = G4bool(*)(const tetra2d&, const tetra2d&);

  G4int FindNeighbor(const std::array<G4int,2>& edge, G4int skip) const;
  G4int FindTetraID(const std::vector<tetra2d>& tetras,
		    const tetra2d& wildTetra, G4int skip,
		    TetraComp tLess) const;
  G4int FirstInteriorTetra();	// Lowest tetra index with all facets shared

  void FindTetrahedron(const G4double point[2], G4double bary[3],
		       G4bool quiet=false) const;

  G4bool Cart2Bary(const G4double point[2], G4double bary[3]) const;
  G4bool BuildT3x2(size_t itet, mat3x2& ET) const;

  G4bool MatInv(const mat2x2& matrix, mat2x2& result, G4bool quiet=false) const;
  G4double BaryNorm(G4double bary[3]) const;
  G4double Det2(const mat2x2& matrix) const;

  // Dump tetrahedron information (neighbors and vertices)
  void PrintTetra(std::ostream& os, G4int iTetra) const;
};

#endif	/* G4CMPBiLinearInterp */
