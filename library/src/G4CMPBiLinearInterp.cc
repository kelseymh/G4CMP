/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20190404  Adapted from BiLinearInterp for use with 2D triangular mesh
// 20190508  Move some 2D/3D common features to new base class
// 20190630  Have MatInv() return error (false), catch up calling chain.
// 20190921  Improve debugging messages, verbose error reports.
// 20190923  Add constructor with neighbors table, use with Clone().
// 20200907  Add FillTInverse() function to precompute invT for Cart2Bary().
//             Use in Cart2Bary() and BuildT4x3() to reduce tracking time.
// 20200908  In MatInv(), clear result first, use new matrix printing.
//		Replace four-arg ctor and UseMesh() with copy constructor.
// 20200914  Include TExtend precalculation in FillTInverse action.
// 20201002  Report tetrahedra errors during FillTInverse() initialization.
// 20240920  G4CMP-244: Replace TetraIdx with function to access G4Cache.

#include "G4CMPBiLinearInterp.hh"
#include "G4CMPConfigManager.hh"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>

using std::array;
using std::map;
using std::sort;
using std::vector;


// Constructors to load mesh from external construction

G4CMPBiLinearInterp::
G4CMPBiLinearInterp(const vector<point2d>& xy, const vector<G4double>& v,
		    const vector<tetra2d>& tetra)
  : G4CMPBiLinearInterp() { UseMesh(xy, v, tetra); }

G4CMPBiLinearInterp::
G4CMPBiLinearInterp(const vector<point3d>& xyz, const vector<G4double>& v,
		    const vector<tetra3d>& tetra)
  : G4CMPBiLinearInterp() { UseMesh(xyz, v, tetra); }

// Copy constructor used by Clone() function

G4CMPBiLinearInterp::G4CMPBiLinearInterp(const G4CMPBiLinearInterp& rhs)
  : G4CMPBiLinearInterp() {
  X = rhs.X;
  V = rhs.V;
  Grad = rhs.Grad;
  Tetrahedra = rhs.Tetrahedra;
  Neighbors = rhs.Neighbors;
  TInverse = rhs.TInverse;
  TInvGood = rhs.TInvGood;
  TExtend  = rhs.TExtend;

  Tetra01 = rhs.Tetra01;	// Not really needed, but for completeness
  Tetra02 = rhs.Tetra02;
  Tetra12 = rhs.Tetra12;

  TetraIdx() = -1;
  TetraStart = rhs.TetraStart;
}


// Load new mesh object and build list of neighbors

void G4CMPBiLinearInterp::UseMesh(const vector<point2d>& xy,
				  const vector<G4double>& v,
				  const vector<tetra2d>& tetra) {
  X = xy;
  V = v;
  Tetrahedra = tetra;
  FillNeighbors();
  FillTInverse();
  FillGradients();
  Initialize();

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
  SaveTetra(savePrefix+"_tetra.dat");
#endif
}


// Load new mesh object using external 3D tables, for client convenience

void G4CMPBiLinearInterp::UseMesh(const vector<point3d>& xyz,
				  const vector<G4double>& v,
				  const vector<tetra3d>& tetra) {
  Compress3DPoints(xyz);
  Compress3DTetras(tetra);
  V = v;
  FillNeighbors();
  FillTInverse();
  FillGradients();
  Initialize();

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
  SaveTetra(savePrefix+"_tetra.dat");
#endif
}


// Return index of tetrahedron with all edges shared, to start FindTetra()

G4int G4CMPBiLinearInterp::FirstInteriorTetra() const {
  G4int minIndex = Neighbors.size()/4;

  for (G4int i=0; i<(G4int)Neighbors.size(); i++) {
    if (*std::min_element(Neighbors[i].begin(), Neighbors[i].end())>minIndex)
      return i;
  }

  return Neighbors.size()/2;
}


// Compress external 3D tables to 2D version (for client convenience)

void G4CMPBiLinearInterp::Compress3DPoints(const std::vector<point3d>& xyz) {
  X.clear();
  X.resize(xyz.size());

  std::transform(xyz.begin(), xyz.end(), X.begin(),
		 [](const point3d& p3d){return point2d{p3d[0],p3d[1]};});
}

void G4CMPBiLinearInterp::Compress3DTetras(const std::vector<tetra3d>& tetra) {
  Tetrahedra.clear();
  Tetrahedra.resize(tetra.size());

  std::transform(tetra.begin(), tetra.end(), Tetrahedra.begin(),
		 [](const tetra3d& t3d){return tetra2d{t3d[0],t3d[1],t3d[2]};});
}


// Tetrahedra sort functions, labelled for each facet option

namespace {
  G4bool tLess01(const tetra2d& a, const tetra2d& b) {
    return ( a[0]<b[0] || (a[0]==b[0] && a[1]<b[1]) );
  }
  
  G4bool tLess02(const tetra2d& a, const tetra2d& b) {
    return ( a[0]<b[0] || (a[0]==b[0] && a[2]<b[2]) );
  }
  
  G4bool tLess12(const tetra2d& a, const tetra2d& b) {
    return ( a[1]<b[1] || (a[1]==b[1] && a[2]<b[2]) );
  }
}

// Process list of defined tetrahedra and build table of neighbors

void G4CMPBiLinearInterp::FillNeighbors() {
  G4cout << "G4CMPBiLinearInterp::FillNeighbors (" << Tetrahedra.size()
	 << " triangles)" << G4endl;

  time_t start, fin;
  std::time(&start);

  // Put the tetrahedra vertices, then the whole list, in indexed order
  for (auto& iTetra: Tetrahedra) sort(iTetra.begin(), iTetra.end());
  sort(Tetrahedra.begin(), Tetrahedra.end());

  // Duplicate list sorted on facets (triplets of vertices)
  Tetra01 = Tetrahedra; sort(Tetra01.begin(), Tetra01.end(), tLess01);
  Tetra02 = Tetrahedra; sort(Tetra02.begin(), Tetra02.end(), tLess02);
  Tetra12 = Tetrahedra; sort(Tetra12.begin(), Tetra12.end(), tLess12);

  G4int Ntet = Tetrahedra.size();		// For convenience below

  Neighbors.clear();
  Neighbors.resize(Ntet, {{-1,-1,-1}});		// Pre-allocate space

  // For each tetrahedron, find another which shares three corners
  for (G4int i=0; i<Ntet; i++) {
    const auto& iTet = Tetrahedra[i];
    Neighbors[i][0] = FindNeighbor({{iTet[1],iTet[2]}}, i);
    Neighbors[i][1] = FindNeighbor({{iTet[0],iTet[2]}}, i);
    Neighbors[i][2] = FindNeighbor({{iTet[0],iTet[1]}}, i);
  }

  std::time(&fin);
  G4cout << "G4CMPBiLinearInterp::FillNeighbors: Took "
         << difftime(fin, start) << " seconds for " << Neighbors.size()
	 << " entries." << G4endl;

}

// Locate other tetrahedron with specified face (excluding "skip" tetrahedron)

G4int G4CMPBiLinearInterp::FindNeighbor(const array<G4int,2>& edge,
					G4int skip) const {
  G4int result = -1;
  result = FindTetraID(Tetra12, {{-1,edge[0],edge[1]}}, skip, tLess12);
  if (result >= 0) return result;	// Successful match

  result = FindTetraID(Tetra02, {{edge[0],-1,edge[1]}}, skip, tLess02);
  if (result >= 0) return result;	// Successful match

  result = FindTetraID(Tetra01, {{edge[0],edge[1],-1}}, skip, tLess01);
  return result;			// If this one failed, they all failed
}

// Locate other tetrahedron with given vertices (excluding "skip" tetrahedron)
// "Wild" means that at least one vertex may be "-1", which matches anything

G4int G4CMPBiLinearInterp::
FindTetraID(const vector<tetra2d>& tetras, const tetra2d& wildTetra, G4int skip,
	    G4CMPBiLinearInterp::TetraComp tLess) const {
  const auto start  = tetras.begin();
  const auto finish = tetras.end();

  // Shared edge (if any) will always be adjacent in sorted table
  auto match = lower_bound(start, finish, wildTetra, tLess);
  if (match == finish) return -1;		// No match at all? PROBLEM!

  G4int index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*match)
		 - Tetrahedra.begin());
  if (index == skip) {				// Move to adjacent entry
    if (match == finish-1) return -1;		// Nothing left to search

    index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*(++match))
	     - Tetrahedra.begin());
  }

  // Test for equality (not < nor >), return index or failure
  return (!tLess(*match,wildTetra) && !tLess(wildTetra,*match)) ? index : -1;
}


// Compute matrices used in tetrahedral barycentric coordinate calculation

void G4CMPBiLinearInterp::FillTInverse() {
#ifdef G4CMPTLI_DEBUG
  G4cout << "G4CMPBiLinearInterp::FillTInverse (" << Tetrahedra.size()
	 << " tetrahedra)" << G4endl;

  time_t start, fin;
  std::time(&start);
#endif

  size_t ntet = Tetrahedra.size();
  TInverse.resize(ntet);		    // Avoid reallocation inside loop
  TExtend.resize(ntet);
  TInvGood.resize(ntet, false);

  mat2x2 T;
  for (size_t itet=0; itet<ntet; itet++) {
    const tetra2d& tetra = Tetrahedra[itet];	// For convenience below
#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 1) {
      G4cout << " Processing Tetrahedra[" << itet << "]: " << tetra << G4endl;
    }
#endif

    for (G4int dim=0; dim<2; ++dim) {
      for (G4int vert=0; vert<2; ++vert) {
	T[dim][vert] = (X[tetra[vert]][dim] - X[tetra[2]][dim]);
      }
    }

    TInvGood[itet] = MatInv(T, TInverse[itet], true);
    BuildT3x2(itet, TExtend[itet]);

    if (!TInvGood[itet]) {
      G4cerr << "ERROR: Non-invertible matrix " << itet << " with " << G4endl;
      for (G4int i=0; i<3; i++) {
	G4cerr << " " << tetra[i] << " @ " << X[tetra[i]] << G4endl;
      }
    }
  }	// for (itet...

#ifdef G4CMPTLI_DEBUG
  std::time(&fin);
  G4cout << "G4CMPBiLinearInterp::FillTInverse: Took "
         << difftime(fin, start) << " seconds for " << TInverse.size()
	 << " entries." << G4endl;
#endif
}


// Compute field (gradient) across each tetrahedron

void G4CMPBiLinearInterp::FillGradients() {
#ifdef G4CMPTLI_DEBUG
  G4cout << "G4CMPBiLinearInterp::FillGradients (" << Tetrahedra.size()
	 << " tetrahedra)" << G4endl;

  time_t start, fin;
  std::time(&start);
#endif

  size_t ntet = Tetrahedra.size();
  Grad.resize(ntet);		    // Avoid reallocation inside loop

  for (size_t itet=0; itet<ntet; itet++) {
    const tetra2d& tetra = Tetrahedra[itet];  // For convenience below
    const mat3x2& ET = TExtend[itet];

    Grad[itet].set((V[tetra[0]]*ET[0][0] + V[tetra[1]]*ET[1][0] +
		    V[tetra[2]]*ET[2][0]),
		   (V[tetra[0]]*ET[0][1] + V[tetra[1]]*ET[1][1] +
		    V[tetra[2]]*ET[2][1]),
		   0.);
#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 1) {
      G4cout << " Computed Grad[" << itet << "]: " << Grad[itet] << G4endl;
    }
#endif
  }	// for (itet...

#ifdef G4CMPTLI_DEBUG
  std::time(&fin);
  G4cout << "G4CMPBiLinearInterp::FillGradients: Took "
         << difftime(fin, start) << " seconds for " << Grad.size()
	 << " entries." << G4endl;
#endif
}


// Evaluate mesh at arbitrary location, returning potential or gradient

G4double 
G4CMPBiLinearInterp::GetValue(const G4double pos[2], G4bool quiet) const {
  G4double bary[3] = { 0. };
  FindTetrahedron(&pos[0], bary, quiet);
    
  if (TetraIdx() == -1) return 0;

  return(V[Tetrahedra[TetraIdx()][0]] * bary[0] +
	 V[Tetrahedra[TetraIdx()][1]] * bary[1] +
	 V[Tetrahedra[TetraIdx()][2]] * bary[2]);    
}

G4ThreeVector 
G4CMPBiLinearInterp::GetGrad(const G4double pos[2], G4bool quiet) const {
  static const G4ThreeVector zero(0.,0.,0.);	// For failure returns

  G4double bary[3] = { 0. };
  FindTetrahedron(pos, bary, quiet);
  return (TetraIdx()<0. ? zero : Grad[TetraIdx()]);
}

void 
G4CMPBiLinearInterp::FindTetrahedron(const G4double pt[2], G4double bary[3],
				      G4bool quiet) const {
  const G4double barySafety = -1e-10;	// Deal with points close to edges

  G4int minBaryIdx = -1;

  G4double bestBary = 0.;	// Norm of barycentric coordinates (below)
  G4int bestTet = -1;

  if (TetraIdx() == -1) TetraIdx() = TetraStart;

#ifdef G4CMPTLI_DEBUG
  if (G4CMPConfigManager::GetVerboseLevel() > 1) {
    G4cout << "FindTetrahedron pt " << pt[0] << " " << pt[1]
	   << "\n starting from TetraIdx " << TetraIdx() << G4endl;
  }
#endif

  // Loop is used to limit search time, does not index tetrahedra
  for (size_t count = 0; count < Tetrahedra.size(); ++count) {
    if (!Cart2Bary(pt,bary)) {	// Get barycentric coord in current tetrahedron
      if (!quiet) {
	G4cerr << "G4CMPBiLinearInterp::FindTetrahedron:"
	       << " Cart2Bary() failed for pt = "
	       << pt[0] << " " << pt[1] << G4endl;

#ifdef G4CMPTLI_DEBUG
	PrintTetra(G4cerr, TetraIdx());
#endif
      }	// if (!quiet)

      TetraIdx() = -1;
      return;
    }	// if (!Cart2bary())

#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 2) {
      G4cout << " Loop " << count << ": Tetra " << TetraIdx() << ": "
	     << Tetrahedra[TetraIdx()] << "\n bary " << bary[0] << " " << bary[1]
	     << " " << bary[2] << " norm " << BaryNorm(bary)
	     << G4endl;
    }
#endif

    // Point is inside current tetrahedron (TetraIdx())
    if (std::all_of(bary, bary+3,
		    [barySafety](G4double b){return b>=barySafety;})) return;

    // Evaluate barycentric distance from current tetrahedron
    G4double newNorm = BaryNorm(bary);
    if (newNorm < bestBary || count == 0) {	// Getting closer
      bestBary = newNorm;
      bestTet  = TetraIdx();

#ifdef G4CMPTLI_DEBUG
      if (G4CMPConfigManager::GetVerboseLevel() > 2) {
	G4cout << " New bestTet " << bestTet << " bestBary = " << bestBary
	       << G4endl;
      }
#endif
    }

    // Point is outside current tetrahedron; shift to nearest neighbor
    minBaryIdx = std::min_element(bary, bary+4) - bary;

    G4int newTetraIdx = Neighbors[TetraIdx()][minBaryIdx];
    if (newTetraIdx == -1) {   // Fell off edge of world
      if (!quiet) {
	G4cerr << "G4CMPBiLinearInterp::FindTetrahedron:"
	       << " Point outside of hull!\n pt = "
	       << pt[0] << " " << pt[1] << G4endl;

#ifdef G4CMPTLI_DEBUG
	PrintTetra(G4cerr, TetraIdx());
#endif
      }

      TetraIdx() = -1;		// Avoids continuing after this
      return;
    }	// if (newTetraIdx == -1)

    TetraIdx() = newTetraIdx;

#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 2) {
      G4cout << " minBaryIdx " << minBaryIdx << ": moved to neighbor tetra "
	     << TetraIdx() << G4endl;
    }
#endif

  }	// for (size_t count=0 ...

  TetraIdx() = bestTet;

#ifdef G4CMPTLI_DEBUG
  if (G4CMPConfigManager::GetVerboseLevel() > 1) {
    Cart2Bary(pt,bary);
    G4cout << "Tetrahedron not found! Using bestTet " << bestTet << " bary "
	   << bary[0] << " " << bary[1] << " " << bary[2]
	   << G4endl;
  }
#endif
}

G4bool 
G4CMPBiLinearInterp::Cart2Bary(const G4double pt[2], G4double bary[3]) const {
  const tetra2d& tetra = Tetrahedra[TetraIdx()]; // For convenience below
  const mat2x2& invT = TInverse[TetraIdx()];
  
  if (TInvGood[TetraIdx()]) {
    bary[2] = 1.0;
    for(G4int k=0; k<2; ++k) {
      bary[k] = (invT[k][0]*(pt[0] - X[tetra[2]][0]) +
		 invT[k][1]*(pt[1] - X[tetra[2]][1]) );
      bary[2] -= bary[k];
    }
  }

  return TInvGood[TetraIdx()];
}

G4double G4CMPBiLinearInterp::BaryNorm(G4double bary[3]) const {
  return (bary[0]*bary[0]+bary[1]*bary[1]+bary[2]*bary[2]);
}

G4bool G4CMPBiLinearInterp::BuildT3x2(size_t itet, mat3x2& ET) const {
  // NOTE:  If matrix inversion failed, invT is set to all zeros
  const mat2x2& invT = TInverse[itet];   // For convenience below
  for (G4int i=0; i<2; ++i) {
    for (G4int j=0; j<2; ++j) {
      ET[i][j] = invT[i][j];
    }
    ET[2][i] = -invT[0][i] - invT[1][i];
  }

  return TInvGood[itet];
}

G4double G4CMPBiLinearInterp::Det2(const mat2x2& matrix) const {
  return (matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]);
}

G4bool G4CMPBiLinearInterp::MatInv(const mat2x2& matrix, mat2x2& result,
				   G4bool quiet) const {
  // Clear any previous result
  std::for_each(result.begin(), result.end(), [](point2d& r){r.fill(0.);});

  // Get determinant and report failure
  G4double determ = Det2(matrix);
  if (!(determ == determ) || fabs(determ) < 1e-9) {
    if (!quiet) {
      G4cerr << "WARNING MatInv got determ " << determ << " zero result"
	     << G4endl;

#ifdef G4CMPTLI_DEBUG
      if (G4CMPConfigManager::GetVerboseLevel() > 1) G4cerr << matrix;
#endif
    }

    return false;
  }

  result[0][0] =  matrix[1][1]/determ;
  result[0][1] = -matrix[0][1]/determ;

  result[1][0] = -matrix[1][0]/determ;
  result[1][1] =  matrix[0][0]/determ;

  return true;
}


// Write data blocks to output file: points and values

void G4CMPBiLinearInterp::SavePoints(const G4String& fname) const {
  G4cout << "Writing points and values to " << fname << G4endl;
  std::ofstream save(fname);
  for (size_t i=0; i<X.size(); i++) {
    save << X[i] << " " << V[i]
	 << std::endl;
  }
}

void G4CMPBiLinearInterp::SaveTetra(const G4String& fname) const {
  G4cout << "Writing tetrahedra and neighbors to " << fname << G4endl;
  std::ofstream save(fname);
    for (size_t i=0; i<Tetrahedra.size(); i++) {
      save << Tetrahedra[i] << "        " << Neighbors[i] << std::endl;
  }
}


// Print out tetrahedral information with coordinates

void G4CMPBiLinearInterp::PrintTetra(std::ostream& os, G4int iTetra) const {
  os << " from tetra " << iTetra << " neighbors " << Neighbors[iTetra] << ":"
     << "\n " << Tetrahedra[iTetra][0] << ": " << X[Tetrahedra[iTetra][0]]
     << "\n " << Tetrahedra[iTetra][1] << ": " << X[Tetrahedra[iTetra][1]]
     << "\n " << Tetrahedra[iTetra][2] << ": " << X[Tetrahedra[iTetra][2]]
     << G4endl;
}
