/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20180525  Protect against "outside of hull" by testing TetraIdx returns;
//	provide "quiet" flag to suppress "outside of hull" messages.
// 20180904  Add constructor to directly load mesh definitions

#include "G4CMPTriLinearInterp.hh"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <algorithm>
#include <ctime>
#include <iostream>

using namespace orgQhull;
using std::array;
using std::map;
using std::sort;
using std::vector;


// Constructors to load mesh and possibly re-triangulate

G4CMPTriLinearInterp::G4CMPTriLinearInterp(const vector<point >& xyz,
					   const vector<G4double>& v)
  : G4CMPTriLinearInterp() { UseMesh(xyz, v); }

G4CMPTriLinearInterp::
G4CMPTriLinearInterp(const vector<point >& xyz, const vector<G4double>& v,
		     const vector<array<G4int,4> >& tetra)
  : G4CMPTriLinearInterp() { UseMesh(xyz, v, tetra); }



// Load new mesh object and possibly re-triangulate

void G4CMPTriLinearInterp::UseMesh(const vector<point > &xyz,
				   const vector<G4double>& v) {
  staleCache = true;
  X = xyz;
  V = v;
  BuildTetraMesh();
  TetraIdx = 0;
}

void G4CMPTriLinearInterp::UseMesh(const vector<point>& xyz,
				   const vector<G4double>& v,
				   const vector<array<G4int,4> >& tetra) {
  staleCache = true;
  X = xyz;
  V = v;
  Tetrahedra = tetra;
  FillNeighbors();
  TetraIdx = 0;
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
  G4bool tLess012(const array<G4int,4>& a, const array<G4int,4>& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[1]<b[1] ||
			     (a[1]==b[1] && a[2]<b[2]))) );
  }
  
  G4bool tLess013(const array<G4int,4>& a, const array<G4int,4>& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[1]<b[1] ||
			     (a[1]==b[1] && a[3]<b[3]))) );
  }
  
  G4bool tLess023(const array<G4int,4>& a, const array<G4int,4>& b) {
    return ( a[0]<b[0] ||
	     (a[0]==b[0] && (a[2]<b[2] ||
			     (a[2]==b[2] && a[3]<b[3]))) );
  }
  
  G4bool tLess123(const array<G4int,4>& a, const array<G4int,4>& b) {
    return ( a[1]<b[1] ||
	     (a[1]==b[1] && (a[2]<b[2] ||
			     (a[2]==b[2] && a[3]<b[3]))) );
  }
}

// Process list of defined tetrahedra and build table of neighbors

void G4CMPTriLinearInterp::FillNeighbors() {
  G4cout << "G4CMPTriLinearInterp::FillNeighbors" << G4endl;

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
         << difftime(fin, start) << " seconds." << G4endl;

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
FindTetraID(const vector<array<G4int,4> >& tetras,
	    const array<G4int,4>& wildTetra, G4int skip,
	    G4CMPTriLinearInterp::TetraComp tLess) const {
  const auto start  = tetras.begin();
  const auto finish = tetras.end();

  // Shared facet (if any) will always be adjacent in sorted table
  auto match = lower_bound(start, finish, wildTetra, tLess);
  if (match == finish) return -1;		// No match at all? PROBLEM!

  G4int index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*match)
		 - Tetrahedra.begin());
  if (index == skip) {				// Move to adjacent entry
    index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*(++match))
	     - Tetrahedra.begin());
  }

  // Test for equality (not < nor >), return index or failure
  return (!tLess(*match,wildTetra) && !tLess(wildTetra,*match)) ? index : -1;
}


// Evaluate mesh at arbitrary location, returning potential or gradient

G4double 
G4CMPTriLinearInterp::GetValue(const G4double pos[3], G4bool quiet) const {
  G4double bary[4] = { 0. };
  FindTetrahedron(&pos[0], bary, quiet);
  staleCache = true;
    
  if (TetraIdx == -1)
    return 0;
  else
    return(V[Tetrahedra[TetraIdx][0]] * bary[0] +
           V[Tetrahedra[TetraIdx][1]] * bary[1] +
           V[Tetrahedra[TetraIdx][2]] * bary[2] +
           V[Tetrahedra[TetraIdx][3]] * bary[3]);    
}

G4ThreeVector 
G4CMPTriLinearInterp::GetGrad(const G4double pos[3], G4bool quiet) const {
  G4double bary[4] = { 0. };
  G4int oldIdx = TetraIdx;
  FindTetrahedron(pos, bary, quiet);

  if (TetraIdx == -1) {
    for (size_t i = 0; i < 3; ++i)
      cachedGrad[i] = 0;
  } else if (TetraIdx != oldIdx || staleCache) {
    G4double ET[4][3];
    BuildT4x3(ET);
    for (size_t i = 0; i < 3; ++i) {
      cachedGrad[i] = V[Tetrahedra[TetraIdx][0]]*ET[0][i] +
                      V[Tetrahedra[TetraIdx][1]]*ET[1][i] +
                      V[Tetrahedra[TetraIdx][2]]*ET[2][i] +
                      V[Tetrahedra[TetraIdx][3]]*ET[3][i];
    }
    staleCache = false;
  }
  return cachedGrad;
}

void 
G4CMPTriLinearInterp::FindTetrahedron(const G4double pt[4], G4double bary[4],
				      G4bool quiet) const {
  const G4double maxError = -1e-10;
  G4int minBaryIdx;
  G4double bestBary[4] = { 0. };
  G4int bestTet = -1;
  if (TetraIdx == -1) TetraIdx = 0;
  for (size_t count = 0; count < Tetrahedra.size(); ++count) {
    Cart2Bary(pt,bary);

    if (bary[3] >= maxError && bary[2] >= maxError &&
	bary[1] >= maxError && bary[0] >= maxError) //bary[3] more likely to be bad.
      return;
    else if ( (bary[0]*bary[0] + bary[1]*bary[1] +
	       bary[2]*bary[2] + bary[3]*bary[3]) < 
	      (bestBary[0]*bestBary[0] + bestBary[1]*bestBary[1] +
	       bestBary[2]*bestBary[2] + bestBary[3]*bestBary[3])
	      || count == 0) {
      for (G4int i = 0; i < 4; ++i)
        bestBary[i] = bary[i];

      bestTet = TetraIdx;
    }

    minBaryIdx = 0;
    for (G4int i = 1; i < 4; ++i)
      if (bary[i] < bary[minBaryIdx]) minBaryIdx = i;

    TetraIdx = Neighbors[TetraIdx][minBaryIdx];
    if (TetraIdx == -1) {
      if (!quiet) {
	G4cerr << "G4CMPTriLinearInterp::FindTetrahedron: Point outside of hull! Check your results."
	       << "\n pt[0] = " << pt[0] << " pt[1] = " << pt[1]
	       << " pt[2] = " << pt[2] << G4endl;
      }

      return;
    }
  }	// for (size_t count=0 ...

  TetraIdx = bestTet;
  Cart2Bary(pt,bary);
  /*G4cout << "G4CMPTriLinearInterp::FindTetrahedron: "
         << "Tetrahedron not found! Using best tetrahedron." << G4endl;
  G4cout << "pt[0] = " << pt[0] << "; pt[1] = " << pt[1] 
         << "; pt[2] = " << pt[2] << ";" << G4endl;
  G4cout << "Best Tetrahedron's barycentric coordinates: " << G4endl;
  G4cout << "bary[0] = " << bary[0] << "; bary[1] = " << bary[1] 
         << "; bary[2] = " << bary[2] << "; bary[3] = " << bary[3] << ";" 
         << G4endl;
         */
}

void G4CMPTriLinearInterp::Cart2Bary(const G4double pt[4], G4double bary[4]) const {
  G4double T[3][3];
  G4double invT[3][3];

  for(G4int dim=0; dim<3; ++dim)
    for(G4int vert=0; vert<3; ++vert)
      T[dim][vert] = X[Tetrahedra[TetraIdx][vert]][dim] - X[Tetrahedra[TetraIdx][3]][dim];

  MatInv(T, invT);

  for(G4int k=0; k<3; ++k)
    bary[k] = invT[k][0]*(pt[0] - X[Tetrahedra[TetraIdx][3]][0]) +
              invT[k][1]*(pt[1] - X[Tetrahedra[TetraIdx][3]][1]) +
              invT[k][2]*(pt[2] - X[Tetrahedra[TetraIdx][3]][2]);

  bary[3] = 1.0 - bary[0] - bary[1] - bary[2];
}

void G4CMPTriLinearInterp::BuildT4x3(G4double ET[4][3]) const {
  G4double T[3][3];
  G4double Tinv[3][3];
  for(G4int dim=0; dim<3; ++dim)
    for(G4int vert=0; vert<3; ++vert)
      T[dim][vert] = X[Tetrahedra[TetraIdx][vert]][dim] - X[Tetrahedra[TetraIdx][3]][dim];

  MatInv(T, Tinv);
  for(G4int i = 0; i < 3; ++i) {
    for(G4int j = 0; j < 3; ++j)
      ET[i][j] = Tinv[i][j];
    ET[3][i] = -Tinv[0][i] - Tinv[1][i] - Tinv[2][i];
  }
}

G4double G4CMPTriLinearInterp::Det3(const G4double matrix[3][3]) const {
  return(matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])
        -matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2])
        +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]));
}

void G4CMPTriLinearInterp::MatInv(const G4double matrix[3][3], G4double result[3][3]) const {
  G4double determ = Det3(matrix);
  result[0][0] = (matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])/determ;
  result[1][0] = (matrix[1][2]*matrix[2][0] - matrix[1][0]*matrix[2][2])/determ;
  result[2][0] = (matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0])/determ;

  result[0][1] = (matrix[0][2]*matrix[2][1] - matrix[0][1]*matrix[2][2])/determ;
  result[1][1] = (matrix[0][0]*matrix[2][2] - matrix[2][0]*matrix[0][2])/determ;
  result[2][1] = (matrix[0][1]*matrix[2][0] - matrix[0][0]*matrix[2][1])/determ;

  result[0][2] = (matrix[0][1]*matrix[1][2] - matrix[1][1]*matrix[0][2])/determ;
  result[1][2] = (matrix[1][0]*matrix[0][2] - matrix[0][0]*matrix[1][2])/determ;
  result[2][2] = (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1])/determ;
}
