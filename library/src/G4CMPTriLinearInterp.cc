/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20180525  Protect against "outside of hull" by testing TetraIdx returns;
//	provide "quiet" flag to suppress "outside of hull" messages.
// 20180904  Add constructor to directly load mesh definitions
// 20180925  Protect memory corruption in passing point coordinates.
// 20180926  Add diagnostic output for debugging field problems.
//		Add starting index for tetrahedral traversal
// 20190226  Provide accessor to replace potentials at mesh points
// 20190404  Change "point" to "point3d" to make way for 2D interpolator.
// 20190508  Move some 2D/3D common features to new base class
// 20190630  Have MatInv() return error (false), catch up calling chain.
// 20190921  Improve debugging messages, verbose error reports.
// 20190923  Add constructor with neighbors table, use with Clone().
// 20200907  Add BuildTInverse() function to precompute invT for Cart2Bary().
//		Use in Cart2Bary() and BuildT4x3() to reduce tracking time.
// 20200908  In MatInv(), clear result first, use new matrix printing.
//		Replace four-arg ctor and UseMesh() with copy constructor.
// 20200914  Include TExtend precalculation in FillTInverse action,
//		gradient (field) precalc in UseMesh functions.
// 20201002  Report tetrahedra errors during FillTInverse() initialization.

#include "G4CMPTriLinearInterp.hh"
#include "G4CMPConfigManager.hh"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>

using namespace orgQhull;
using std::array;
using std::map;
using std::sort;
using std::vector;


// Constructors to load mesh and possibly re-triangulate

G4CMPTriLinearInterp::G4CMPTriLinearInterp(const vector<point3d>& xyz,
					   const vector<G4double>& v)
  : G4CMPTriLinearInterp() { UseMesh(xyz, v); }

G4CMPTriLinearInterp::
G4CMPTriLinearInterp(const vector<point3d>& xyz, const vector<G4double>& v,
		     const vector<tetra3d>& tetra)
  : G4CMPTriLinearInterp() { UseMesh(xyz, v, tetra); }

// Copy constructor used by Clone() function

G4CMPTriLinearInterp::G4CMPTriLinearInterp(const G4CMPTriLinearInterp& rhs)
  : G4CMPTriLinearInterp() {
  X = rhs.X;
  V = rhs.V;
  Grad = rhs.Grad;
  Tetrahedra = rhs.Tetrahedra;
  Neighbors = rhs.Neighbors;
  TInverse = rhs.TInverse;
  TInvGood = rhs.TInvGood;
  TExtend  = rhs.TExtend;

  Tetra012 = rhs.Tetra012;	// Not really needed, but for completeness
  Tetra013 = rhs.Tetra013;
  Tetra023 = rhs.Tetra023;
  Tetra123 = rhs.Tetra123;

  TetraIdx = -1;
  TetraStart = rhs.TetraStart;
}


// Load new mesh object and possibly re-triangulate

void G4CMPTriLinearInterp::UseMesh(const vector<point3d> &xyz,
				   const vector<G4double>& v) {
  X = xyz;
  V = v;
  BuildTetraMesh();
  FillTInverse();
  FillGradients();

  TetraIdx = -1;
  TetraStart = FirstInteriorTetra();

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
  SaveTetra(savePrefix+"_tetra.dat");
#endif
}

void G4CMPTriLinearInterp::UseMesh(const vector<point3d>& xyz,
				   const vector<G4double>& v,
				   const vector<tetra3d>& tetra) {
  X = xyz;
  V = v;
  Tetrahedra = tetra;
  FillNeighbors();
  FillTInverse();
  FillGradients();

  TetraIdx = -1;
  TetraStart = FirstInteriorTetra();

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
  SaveTetra(savePrefix+"_tetra.dat");
#endif
}


// Generate new Delaunay triagulation for current mesh of points

void G4CMPTriLinearInterp::BuildTetraMesh() {
  time_t start, fin;
  G4cout << "G4CMPTriLinearInterp::Constructor: Creating Tetrahedral Mesh..."
         << G4endl;
  std::time(&start);
  /* Qhull requires a column-major array of the
   * 3D points. i.e., [x1,y1,z1,x2,y2,z2,...]
   */
  G4double* boxPoints = new G4double[3*X.size()];
      
  for (size_t i=0, e=X.size(); i<e; ++i) {
    boxPoints[i*3] = X[i][0];
    boxPoints[i*3+1]= X[i][1];
    boxPoints[i*3+2]= X[i][2];
  }
    
  /* Run Qhull
   * First parameter is a command to generate random points. We don't want that.
   * Second parameter is dimension of points (3)
   * Third parameter is number of points
   * Forth parameter is the 1D column-major array of points
   * Fifth parameter is qhull option flags.
   *   d = delaunay triangulation
   *   Qt = triangulate output (makes all regions simplicial)
   *   Qbb = Scales the paraboloid that Qhull creates. This helps with precision
   *   Qz = Add a point at infinity. This somehow helps with precision...
   */
  Qhull hull = Qhull("", 3, X.size(), boxPoints, "d Qt Qz Qbb");
        
  QhullFacet facet, neighbor;
  QhullVertex vertex;
  QhullLinkedList<QhullFacet>::iterator fItr;
  QhullSet<QhullFacet>::iterator nItr;
  QhullSet<QhullVertex>::iterator vItr;
  map<G4int, G4int> ID2Idx;
  G4int numTet = 0, j;
  vector<array<G4int, 4> > tmpTetrahedra(hull.facetCount(), {{0,0,0,0}});
  vector<array<G4int, 4> > tmpNeighbors(hull.facetCount(), {{-1,-1,-1,-1}});
  for (fItr = hull.facetList().begin();fItr != hull.facetList().end(); fItr++) {
    facet = *fItr;
    if (!facet.isUpperDelaunay()) {
      ID2Idx[facet.id()] = numTet;
      j = 0;
      for (vItr = facet.vertices().begin(); vItr != facet.vertices().end(); vItr++) {
        vertex = *vItr;
        tmpTetrahedra[numTet][j++] = FindPointID(vertex.point().toStdVector(), vertex.id());
      }
      j = 0;
      for (nItr = facet.neighborFacets().begin(); nItr != facet.neighborFacets().end(); nItr++) {
        neighbor = *nItr;
        if (!neighbor.isUpperDelaunay())
          tmpNeighbors[numTet][j] = neighbor.id();
        ++j;
      }
      ++numTet;
    }
  }
    
  for (G4int i = 0; i < numTet; ++i)
    for (j = 0; j < 4; ++j) {
      if (tmpNeighbors[i][j] != -1)
        tmpNeighbors[i][j] = ID2Idx[tmpNeighbors[i][j]];
    }

  tmpTetrahedra.resize(numTet);
  tmpNeighbors.resize(numTet);

  Tetrahedra.swap(tmpTetrahedra);
  Neighbors.swap(tmpNeighbors);

  delete[] boxPoints;

  std::time(&fin);
  G4cout << "G4CMPTriLinearInterp::Constructor: Took "
         << difftime(fin, start) << " seconds." << G4endl;
}

// Get index of specified mesh point (no interpolation!)

G4int G4CMPTriLinearInterp::FindPointID(const vector<G4double>& pt,
                                        const G4int id) const {
  if (qhull2x.count(id)) {
    return qhull2x[id];
  }
	
  G4int min = 0, max = X.size()-1, mid = X.size()/2;
  G4int d = 0;
  while(true) {
    if (X[mid][0] == pt[0]) {
      d = 0;
      while(X[mid+d][0] == pt[0] || X[mid-d][0] == pt[0]) {
        if (X[mid+d][0] == pt[0] && X[mid+d][1] == pt[1] && X[mid+d][2] == pt[2]) {
          qhull2x[id] = mid+d;
          return mid + d;
        } else if (X[mid-d][0] == pt[0] && X[mid-d][1] == pt[1] && X[mid-d][2] == pt[2]) {
          qhull2x[id] = mid-d;
          return mid - d;
        }
        ++d;
      }
      while(X[mid-d][0] == pt[0]) {
        if (X[mid-d][0] == pt[0] && X[mid-d][1] == pt[1] && X[mid-d][2] == pt[2]) {
          qhull2x[id] = mid-d;
          return mid - d;
        }
        ++d;
      }
    } else if (X[mid][0] < pt[0]) {
      min = mid + 1;
      mid = min + (G4int)ceil((max - min)/2.);
    } else {
      max = mid - 1;
      mid = min + (G4int)floor((max - min)/2.);
    }
  }
}


// Tetrahedra sort functions, labelled for each facet option

namespace {
  G4bool tLess012(const tetra3d& a, const tetra3d& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[1]<b[1] ||
			     (a[1]==b[1] && a[2]<b[2]))) );
  }
  
  G4bool tLess013(const tetra3d& a, const tetra3d& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[1]<b[1] ||
			     (a[1]==b[1] && a[3]<b[3]))) );
  }
  
  G4bool tLess023(const tetra3d& a, const tetra3d& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[2]<b[2] ||
			     (a[2]==b[2] && a[3]<b[3]))) );
  }
  
  G4bool tLess123(const tetra3d& a, const tetra3d& b) {
    return ( a[1]<b[1] ||
	     (a[1]==b[1] && (a[2]<b[2] ||
			     (a[2]==b[2] && a[3]<b[3]))) );
  }
}

// Process list of defined tetrahedra and build table of neighbors

void G4CMPTriLinearInterp::FillNeighbors() {
  G4cout << "G4CMPTriLinearInterp::FillNeighbors (" << Tetrahedra.size()
	 << " tetrahedra)" << G4endl;

  time_t start, fin;
  std::time(&start);

  // Put the tetrahedra vertices, then the whole list, in indexed order
  for (auto& iTetra: Tetrahedra) sort(iTetra.begin(), iTetra.end());
  sort(Tetrahedra.begin(), Tetrahedra.end());

  // Duplicate list sorted on facets (triplets of vertices)
  Tetra012 = Tetrahedra; sort(Tetra012.begin(), Tetra012.end(), tLess012);
  Tetra013 = Tetrahedra; sort(Tetra013.begin(), Tetra013.end(), tLess013);
  Tetra023 = Tetrahedra; sort(Tetra023.begin(), Tetra023.end(), tLess023);
  Tetra123 = Tetrahedra; sort(Tetra123.begin(), Tetra123.end(), tLess123);

  G4int Ntet = Tetrahedra.size();		// For convenience below

  Neighbors.clear();
  Neighbors.resize(Ntet, {{-1,-1,-1,-1}});	// Pre-allocate space

  // For each tetrahedron, find another which shares three corners
  for (G4int i=0; i<Ntet; i++) {
    const auto& iTet = Tetrahedra[i];
    Neighbors[i][0] = FindNeighbor({{iTet[1],iTet[2],iTet[3]}}, i);
    Neighbors[i][1] = FindNeighbor({{iTet[0],iTet[2],iTet[3]}}, i);
    Neighbors[i][2] = FindNeighbor({{iTet[0],iTet[1],iTet[3]}}, i);
    Neighbors[i][3] = FindNeighbor({{iTet[0],iTet[1],iTet[2]}}, i);
  }

  std::time(&fin);
  G4cout << "G4CMPTriLinearInterp::FillNeighbors: Took "
         << difftime(fin, start) << " seconds for " << Neighbors.size()
	 << " entries." << G4endl;

}

// Locate other tetrahedron with specified face (excluding "skip" tetrahedron)

G4int G4CMPTriLinearInterp::FindNeighbor(const array<G4int,3>& facet,
					G4int skip) const {
  G4int result = -1;
  result = FindTetraID(Tetra123, {{-1,facet[0],facet[1],facet[2]}}, skip, tLess123);
  if (result >= 0) return result;	// Successful match

  result = FindTetraID(Tetra023, {{facet[0],-1,facet[1],facet[2]}}, skip, tLess023);
  if (result >= 0) return result;	// Successful match

  result = FindTetraID(Tetra013, {{facet[0],facet[1],-1,facet[2]}}, skip, tLess013);
  if (result >= 0) return result;	// Successful match

  result = FindTetraID(Tetra012, {{facet[0],facet[1],facet[2],-1}}, skip, tLess012);
  return result;			// If this one failed, they all failed
}

// Locate other tetrahedron with given vertices (excluding "skip" tetrahedron)
// "Wild" means that at least one vertex may be "-1", which matches anything

G4int G4CMPTriLinearInterp::
FindTetraID(const vector<tetra3d>& tetras, const tetra3d& wildTetra, G4int skip,
	    G4CMPTriLinearInterp::TetraComp tLess) const {
  const auto start  = tetras.begin();
  const auto finish = tetras.end();

  // Shared facet (if any) will always be adjacent in sorted table
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

void G4CMPTriLinearInterp::FillTInverse() {
#ifdef G4CMPTLI_DEBUG
  G4cout << "G4CMPTriLinearInterp::FillTInverse (" << Tetrahedra.size()
	 << " tetrahedra)" << G4endl;

  time_t start, fin;
  std::time(&start);
#endif

  size_t ntet = Tetrahedra.size();
  TInverse.resize(ntet);		    // Avoid reallocation inside loop
  TExtend.resize(ntet);
  TInvGood.resize(ntet, false);

  mat3x3 T;
  for (size_t itet=0; itet<ntet; itet++) {
    const tetra3d& tetra = Tetrahedra[itet];	// For convenience below
#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 1) {
      G4cout << " Processing Tetrahedra[" << itet << "]: " << tetra << G4endl;
    }
#endif

    for (G4int dim=0; dim<3; ++dim) {
      for (G4int vert=0; vert<3; ++vert) {
	T[dim][vert] = (X[tetra[vert]][dim] - X[tetra[3]][dim]);
      }
    }

    TInvGood[itet] = MatInv(T, TInverse[itet], true);
    BuildT4x3(itet, TExtend[itet]);

    if (!TInvGood[itet]) {
      G4cerr << "ERROR: Non-invertible matrix " << itet << " with " << G4endl;
      for (G4int i=0; i<4; i++) {
	G4cerr << " " << tetra[i] << " @ " << X[tetra[i]] << G4endl;
      }
    }
  }	// for (itet...

#ifdef G4CMPTLI_DEBUG
  std::time(&fin);
  G4cout << "G4CMPTriLinearInterp::FillTInverse: Took "
         << difftime(fin, start) << " seconds for " << TInverse.size()
	 << " entries." << G4endl;
#endif
}


// Compute field (gradient) across each tetrahedron

void G4CMPTriLinearInterp::FillGradients() {
#ifdef G4CMPTLI_DEBUG
  G4cout << "G4CMPTriLinearInterp::FillGradients (" << Tetrahedra.size()
	 << " tetrahedra)" << G4endl;

  time_t start, fin;
  std::time(&start);
#endif

  size_t ntet = Tetrahedra.size();
  Grad.resize(ntet);		    // Avoid reallocation inside loop

  for (size_t itet=0; itet<ntet; itet++) {
    const tetra3d& tetra = Tetrahedra[itet];  // For convenience below
    const mat4x3& ET = TExtend[itet];

    Grad[itet].set((V[tetra[0]]*ET[0][0] + V[tetra[1]]*ET[1][0] +
		    V[tetra[2]]*ET[2][0] + V[tetra[3]]*ET[3][0]),
		   (V[tetra[0]]*ET[0][1] + V[tetra[1]]*ET[1][1] +
		    V[tetra[2]]*ET[2][1] + V[tetra[3]]*ET[3][1]),
		   (V[tetra[0]]*ET[0][2] + V[tetra[1]]*ET[1][2] +
		    V[tetra[2]]*ET[2][2] + V[tetra[3]]*ET[3][2])
		   );
#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 1) {
      G4cout << " Computed Grad[" << itet << "]: " << Grad[itet] << G4endl;
    }
#endif
  }	// for (itet...

#ifdef G4CMPTLI_DEBUG
  std::time(&fin);
  G4cout << "G4CMPTriLinearInterp::FillGradients: Took "
         << difftime(fin, start) << " seconds for " << Grad.size()
	 << " entries." << G4endl;
#endif
}


// Return index of tetrahedron with all facets shared, to start FindTetra()

G4int G4CMPTriLinearInterp::FirstInteriorTetra() {
  G4int minIndex = Neighbors.size()/4;

  for (G4int i=0; i<(G4int)Neighbors.size(); i++) {
    if (*std::min_element(Neighbors[i].begin(), Neighbors[i].end())>minIndex)
      return i;
  }

  return Neighbors.size()/2;
}

// Evaluate mesh at arbitrary location, returning potential or gradient

G4double 
G4CMPTriLinearInterp::GetValue(const G4double pos[3], G4bool quiet) const {
  G4double bary[4] = { 0. };
  FindTetrahedron(&pos[0], bary, quiet);
    
  if (TetraIdx == -1) return 0;

  return(V[Tetrahedra[TetraIdx][0]] * bary[0] +
	 V[Tetrahedra[TetraIdx][1]] * bary[1] +
	 V[Tetrahedra[TetraIdx][2]] * bary[2] +
	 V[Tetrahedra[TetraIdx][3]] * bary[3]);    
}

G4ThreeVector 
G4CMPTriLinearInterp::GetGrad(const G4double pos[3], G4bool quiet) const {
  static const G4ThreeVector zero(0.,0.,0.);	// For failure returns

  G4double bary[4] = { 0. };
  FindTetrahedron(pos, bary, quiet);
  return (TetraIdx<0. ? zero : Grad[TetraIdx]);
}


// Identify tetrahedron enclosing point, returning barycentric coords

void 
G4CMPTriLinearInterp::FindTetrahedron(const G4double pt[3], G4double bary[4],
				      G4bool quiet) const {
  const G4double barySafety = -1e-10;	// Deal with points close to facets

  G4double bestBary = 0.;	// Norm of barycentric coordinates (below)
  G4int bestTet = -1;

  if (TetraIdx == -1) TetraIdx = TetraStart;

#ifdef G4CMPTLI_DEBUG
  if (G4CMPConfigManager::GetVerboseLevel() > 1) {
    G4cout << "FindTetrahedron pt " << pt[0] << " " << pt[1] << " " << pt[2]
	   << "\n starting from TetraIdx " << TetraIdx << G4endl;
  }
#endif

  // Loop is used to limit search time, does not index tetrahedra
  for (size_t count = 0; count < Tetrahedra.size(); ++count) {
    if (!Cart2Bary(pt,bary)) {	// Get barycentric coord in current tetrahedron
      if (!quiet) {
	G4cerr << "G4CMPTriLinearInterp::FindTetrahedron:"
	       << " Cart2Bary() failed for pt = "
	       << pt[0] << " " << pt[1] << " " << pt[2] << G4endl;

#ifdef G4CMPTLI_DEBUG
	PrintTetra(G4cerr, TetraIdx);
#endif
      }	// if (!quiet)

      TetraIdx = -1;
      return;
    }	// if (!Cart2bary())

#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 2) {
      G4cout << " Loop " << count << ": Tetra " << TetraIdx << ": "
	     << Tetrahedra[TetraIdx] << "\n bary " << bary[0] << " " << bary[1]
	     << " " << bary[2] << " " << bary[3] << " norm " << BaryNorm(bary)
	     << G4endl;
    }
#endif

    // Point is inside current tetrahedron (TetraIdx)
    if (std::all_of(bary, bary+4,
		    [barySafety](G4double b){return b>=barySafety;})) return;

    // Evaluate barycentric distance from current tetrahedron
    G4double newNorm = BaryNorm(bary);
    if (newNorm < bestBary || count == 0) {	// Getting closer
      bestBary = newNorm;
      bestTet  = TetraIdx;

#ifdef G4CMPTLI_DEBUG
      if (G4CMPConfigManager::GetVerboseLevel() > 2) {
	G4cout << " New bestTet " << bestTet << " bestBary = " << bestBary
	       << G4endl;
      }
#endif
    }

    // Point is outside current tetrahedron; shift to nearest neighbor
    G4int minBaryIdx = std::min_element(bary, bary+4) - bary;

    G4int newTetraIdx = Neighbors[TetraIdx][minBaryIdx];
    if (newTetraIdx == -1) {	// Fell off edge of world
      if (!quiet) {
	G4cerr << "G4CMPTriLinearInterp::FindTetrahedron:"
	       << " Point outside of hull!\n pt = "
	       << pt[0] << " " << pt[1] << " " << pt[2] << G4endl;

#ifdef G4CMPTLI_DEBUG
	PrintTetra(G4cerr, TetraIdx);
#endif
      }

      TetraIdx = -1;		// Avoids continuing after this
      return;
    }	// if (newTetraIdx == -1)

    TetraIdx = newTetraIdx;

#ifdef G4CMPTLI_DEBUG
    if (G4CMPConfigManager::GetVerboseLevel() > 2) {
      G4cout << " minBaryIdx " << minBaryIdx << ": moved to neighbor tetra "
	     << TetraIdx << G4endl;
    }
#endif

  }	// for (size_t count=0 ...

  TetraIdx = bestTet;

#ifdef G4CMPTLI_DEBUG
  if (G4CMPConfigManager::GetVerboseLevel() > 1) {
    Cart2Bary(pt,bary);
    G4cout << "Tetrahedron not found! Using bestTet " << bestTet << " bary "
	   << bary[0] << " " << bary[1] << " " << bary[2] << " " << bary[3]
	   << G4endl;
  }
#endif
}

G4bool
G4CMPTriLinearInterp::Cart2Bary(const G4double pt[3], G4double bary[4]) const {
  const tetra3d& tetra = Tetrahedra[TetraIdx];	// For convenience below
  const mat3x3& invT = TInverse[TetraIdx];

  if (TInvGood[TetraIdx]) {
    bary[3] = 1.0;
    for(G4int k=0; k<3; ++k) {
      bary[k] = (invT[k][0]*(pt[0] - X[tetra[3]][0]) +
		 invT[k][1]*(pt[1] - X[tetra[3]][1]) +
		 invT[k][2]*(pt[2] - X[tetra[3]][2]) );
      bary[3] -= bary[k];
    }
  }

  return TInvGood[TetraIdx];
}

G4double G4CMPTriLinearInterp::BaryNorm(G4double bary[4]) const {
  return (bary[0]*bary[0]+bary[1]*bary[1]+bary[2]*bary[2]+bary[3]*bary[3]);
}

G4bool G4CMPTriLinearInterp::BuildT4x3(size_t iTet, mat4x3& ET) const {
  // NOTE:  If matrix inversion failed, invT is set to all zeros
  const mat3x3& invT = TInverse[iTet];	// For convenience below
  for (G4int i=0; i<3; ++i) {
    for (G4int j=0; j<3; ++j) {
      ET[i][j] = invT[i][j];
    }
    ET[3][i] = -invT[0][i] - invT[1][i] - invT[2][i];
  }

  return TInvGood[iTet];
}

G4double G4CMPTriLinearInterp::Det3(const mat3x3& matrix) const {
  return(matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])
        +matrix[0][1]*(matrix[1][2]*matrix[2][0]-matrix[2][2]*matrix[1][0])
        +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]));
}

G4bool G4CMPTriLinearInterp::MatInv(const mat3x3& matrix, mat3x3& result,
				    G4bool quiet) const {
  // Clear any previous result
  std::for_each(result.begin(), result.end(), [](point3d& r){r.fill(0.);});

  // Get determinant and report failure
  G4double determ = Det3(matrix);
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

  result[0][0] = (matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])/determ;
  result[1][0] = (matrix[1][2]*matrix[2][0] - matrix[1][0]*matrix[2][2])/determ;
  result[2][0] = (matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0])/determ;

  result[0][1] = (matrix[0][2]*matrix[2][1] - matrix[0][1]*matrix[2][2])/determ;
  result[1][1] = (matrix[0][0]*matrix[2][2] - matrix[2][0]*matrix[0][2])/determ;
  result[2][1] = (matrix[0][1]*matrix[2][0] - matrix[0][0]*matrix[2][1])/determ;

  result[0][2] = (matrix[0][1]*matrix[1][2] - matrix[1][1]*matrix[0][2])/determ;
  result[1][2] = (matrix[1][0]*matrix[0][2] - matrix[0][0]*matrix[1][2])/determ;
  result[2][2] = (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1])/determ;

  return true;
}


// Write data blocks to output file: points and values

void G4CMPTriLinearInterp::SavePoints(const G4String& fname) const {
  G4cout << "Writing points and values to " << fname << G4endl;
  std::ofstream save(fname);
  for (size_t i=0; i<X.size(); i++) {
    save << X[i] << " " << V[i]
	 << std::endl;
  }
}

void G4CMPTriLinearInterp::SaveTetra(const G4String& fname) const {
  G4cout << "Writing tetrahedra and neighbors to " << fname << G4endl;
  std::ofstream save(fname);
    for (size_t i=0; i<Tetrahedra.size(); i++) {
      save << Tetrahedra[i] << "        " << Neighbors[i] << std::endl;
  }
}


// Print out tetrahedral information with coordinates

void G4CMPTriLinearInterp::PrintTetra(std::ostream& os, G4int iTetra) const {
  os << " from tetra " << iTetra << " neighbors " << Neighbors[iTetra] << ":"
     << "\n " << Tetrahedra[iTetra][0] << ": " << X[Tetrahedra[iTetra][0]]
     << "\n " << Tetrahedra[iTetra][1] << ": " << X[Tetrahedra[iTetra][1]]
     << "\n " << Tetrahedra[iTetra][2] << ": " << X[Tetrahedra[iTetra][2]]
     << "\n " << Tetrahedra[iTetra][3] << ": " << X[Tetrahedra[iTetra][3]]
     << G4endl;
}
