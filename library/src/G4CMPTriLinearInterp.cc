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

#include "G4CMPTriLinearInterp.hh"
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



// Load new mesh object and possibly re-triangulate

void G4CMPTriLinearInterp::UseMesh(const vector<point3d> &xyz,
				   const vector<G4double>& v) {
  staleCache = true;
  X = xyz;
  V = v;
  BuildTetraMesh();
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
  staleCache = true;
  X = xyz;
  V = v;
  Tetrahedra = tetra;
  FillNeighbors();
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
    index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*(++match))
	     - Tetrahedra.begin());
  }

  // Test for equality (not < nor >), return index or failure
  return (!tLess(*match,wildTetra) && !tLess(wildTetra,*match)) ? index : -1;
}


// Return index of tetrahedron with all facets shared, to start FindTetra()

G4int G4CMPTriLinearInterp::FirstInteriorTetra() {
  G4int minIndex = Neighbors.size()/4;

  for (G4int i=0; i<(G4int)Neighbors.size(); i++) {
    if (Neighbors[i][0]>minIndex && Neighbors[i][1]>minIndex &&
	Neighbors[i][2]>minIndex && Neighbors[i][3]>minIndex) return i;
  }

  return Neighbors.size()/2;
}

// Evaluate mesh at arbitrary location, returning potential or gradient

G4double 
G4CMPTriLinearInterp::GetValue(const G4double pos[3], G4bool quiet) const {
  G4double bary[4] = { 0. };
  FindTetrahedron(&pos[0], bary, quiet);
  staleCache = true;
    
  if (TetraIdx == -1) return 0;
  else {
    return(V[Tetrahedra[TetraIdx][0]] * bary[0] +
           V[Tetrahedra[TetraIdx][1]] * bary[1] +
           V[Tetrahedra[TetraIdx][2]] * bary[2] +
           V[Tetrahedra[TetraIdx][3]] * bary[3]);    
  }
}

G4ThreeVector 
G4CMPTriLinearInterp::GetGrad(const G4double pos[3], G4bool quiet) const {
  G4double bary[4] = { 0. };
  G4int oldIdx = TetraIdx;
  FindTetrahedron(pos, bary, quiet);

  if (TetraIdx == -1) cachedGrad.set(0.,0.,0.);
  else if (TetraIdx != oldIdx || staleCache) {
    G4double ET[4][3];
    if (!BuildT4x3(ET)) cachedGrad.set(0.,0.,0.);	// Failed MatInv()
    else {
      for (size_t i = 0; i < 3; ++i) {
	cachedGrad[i] = V[Tetrahedra[TetraIdx][0]]*ET[0][i] +
                        V[Tetrahedra[TetraIdx][1]]*ET[1][i] +
                        V[Tetrahedra[TetraIdx][2]]*ET[2][i] +
                        V[Tetrahedra[TetraIdx][3]]*ET[3][i];
      }
      staleCache = false;
    }
  }

  return cachedGrad;
}

void 
G4CMPTriLinearInterp::FindTetrahedron(const G4double pt[3], G4double bary[4],
				      G4bool quiet) const {
  const G4double barySafety = -1e-10;	// Deal with points close to facets

  G4int minBaryIdx = -1;

  G4double bestBary = 0.;	// Norm of barycentric coordinates (below)
  G4int bestTet = -1;

  if (TetraIdx == -1) TetraIdx = TetraStart;

#ifdef G4CMPTLI_DEBUG
  G4cout << "FindTetrahedron pt " << pt[0] << " " << pt[1] << " " << pt[2]
	 << "\n starting from TetraIdx " << TetraIdx << G4endl;
#endif

  for (size_t count = 0; count < Tetrahedra.size(); ++count) {
    if (!Cart2Bary(pt,bary)) {	// Get barycentric coord in current tetrahedron
      if (!quiet) {
	G4cerr << "G4CMPTriLinearInterp::FindTetrahedron:"
	       << " Cart2Bary() failed for pt = "
	       << pt[0] << " " << pt[1] << " " << pt[2] << G4endl;
      }

      TetraIdx = -1;
      return;
    }

#ifdef G4CMPTLI_DEBUG
    G4cout << " Loop " << count << ": Tetra " << TetraIdx << ": "
	   << Tetrahedra[TetraIdx] << "\n bary " << bary[0] << " " << bary[1]
	   << " " << bary[2] << " " << bary[3] << " norm " << BaryNorm(bary)
	   << G4endl;
#endif

    // Point is inside current tetrahedron (TetraIdx)
    if (std::all_of(bary, bary+4,
		    [barySafety](G4double b){return b>=barySafety;})) return;

    // Evaluate barycentric distance from current tetrahedron
    if (BaryNorm(bary) < bestBary || count == 0) {	// Getting closer
      bestBary = BaryNorm(bary);
      bestTet  = TetraIdx;
#ifdef G4CMPTLI_DEBUG
      G4cout << " New bestTet " << bestTet << " bestBary = " << bestBary
	     << G4endl;
#endif
    }

    // Point is outside current tetrahedron; shift to nearest neighbor
    minBaryIdx = 0;
    for (G4int i=1; i<4; ++i)
      if (bary[i] < bary[minBaryIdx]) minBaryIdx = i;

    TetraIdx = Neighbors[TetraIdx][minBaryIdx];

#ifdef G4CMPTLI_DEBUG
    G4cout << " minBaryIdx " << minBaryIdx << ": moved to neighbor tetra "
	   << TetraIdx << G4endl;
#endif
    if (TetraIdx == -1) {
      if (!quiet) {
	G4cerr << "G4CMPTriLinearInterp::FindTetrahedron:"
	       << " Point outside of hull!\n pt = "
	       << pt[0] << " " << pt[1] << " " << pt[2] << G4endl;
      }
      return;
    }
  }	// for (size_t count=0 ...

  TetraIdx = bestTet;
  Cart2Bary(pt,bary);		// Don't need to check return; succeeded above

#ifdef G4CMPTLI_DEBUG
  G4cout << "Tetrahedron not found! Using bestTet " << bestTet << " bary "
	 << bary[0] << " " << bary[1] << " " << bary[2] << " " << bary[3]
         << G4endl;
#endif
}

G4bool
G4CMPTriLinearInterp::Cart2Bary(const G4double pt[3], G4double bary[4]) const {
  G4double T[3][3];
  G4double invT[3][3];

  for(G4int dim=0; dim<3; ++dim)
    for(G4int vert=0; vert<3; ++vert)
      T[dim][vert] = (X[Tetrahedra[TetraIdx][vert]][dim]
		      - X[Tetrahedra[TetraIdx][3]][dim]);

  G4bool goodInv = MatInv(T, invT);

  for(G4int k=0; k<3; ++k) {
    bary[k] = invT[k][0]*(pt[0] - X[Tetrahedra[TetraIdx][3]][0]) +
              invT[k][1]*(pt[1] - X[Tetrahedra[TetraIdx][3]][1]) +
              invT[k][2]*(pt[2] - X[Tetrahedra[TetraIdx][3]][2]);
  }

  bary[3] = 1.0 - bary[0] - bary[1] - bary[2];

  return goodInv;
}

G4double G4CMPTriLinearInterp::BaryNorm(G4double bary[4]) const {
  return (bary[0]*bary[0]+bary[1]*bary[1]+bary[2]*bary[2]+bary[3]*bary[3]);
}

G4bool G4CMPTriLinearInterp::BuildT4x3(G4double ET[4][3]) const {
  G4double T[3][3];
  G4double Tinv[3][3];
  for (G4int dim=0; dim<3; ++dim)
    for (G4int vert=0; vert<3; ++vert)
      T[dim][vert] = (X[Tetrahedra[TetraIdx][vert]][dim]
		      - X[Tetrahedra[TetraIdx][3]][dim]);

  G4bool goodInv = MatInv(T, Tinv);
  if (goodInv) {
    for (G4int i=0; i<3; ++i) {
      for (G4int j=0; j<3; ++j)
	ET[i][j] = Tinv[i][j];
      ET[3][i] = -Tinv[0][i] - Tinv[1][i] - Tinv[2][i];
    }
  } else {
    for (G4int i=0; i<4; ++i) {
      for (G4int j=0; j<3; ++j) {
	ET[i][j] = 0.;
      }
    }
  }

  return goodInv;
}

G4double G4CMPTriLinearInterp::Det3(const G4double matrix[3][3]) const {
  return(matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])
        -matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2])
        +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]));
}

G4bool G4CMPTriLinearInterp::MatInv(const G4double matrix[3][3],
				    G4double result[3][3]) const {
  G4double determ = Det3(matrix);
  if (!(determ == determ) || fabs(determ) < 1e-9) {
    G4cerr << "WARNING MatInv got determ " << determ << " zero result" << G4endl;
#ifdef G4CMPTLI_DEBUG
    G4cerr << " "  << matrix[0][0] << " " << matrix[0][1] << " " << matrix[0][2]
	   <<"\n " << matrix[1][0] << " " << matrix[1][1] << " " << matrix[1][2]
	   <<"\n " << matrix[2][0] << " " << matrix[2][1] << " " << matrix[2][2]
	   << G4endl;
#endif

    for (size_t i=0; i<3; i++) {
      for (size_t j=0; j<3; j++) {
	result[i][j] = 0.;
      }
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
