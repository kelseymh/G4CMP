/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20190404  Adapted from BiLinearInterp for use with 2D triangular mesh
// 20190508  Move some 2D/3D common features to new base class

#include "G4CMPBiLinearInterp.hh"
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


// Load new mesh object and build list of neighbors

void G4CMPBiLinearInterp::UseMesh(const vector<point2d>& xy,
				  const vector<G4double>& v,
				  const vector<tetra2d>& tetra) {
  staleCache = true;
  X = xy;
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

// Load new mesh object using external 3D tables, for client convenience

void G4CMPBiLinearInterp::UseMesh(const vector<point3d>& xyz,
				  const vector<G4double>& v,
				  const vector<tetra3d>& tetra) {
  staleCache = true;
  Compress3DPoints(xyz);
  Compress3DTetras(tetra);
  V = v;
  FillNeighbors();
  TetraIdx = -1;
  TetraStart = FirstInteriorTetra();

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
  SaveTetra(savePrefix+"_tetra.dat");
#endif
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
    index = (lower_bound(Tetrahedra.begin(),Tetrahedra.end(),*(++match))
	     - Tetrahedra.begin());
  }

  // Test for equality (not < nor >), return index or failure
  return (!tLess(*match,wildTetra) && !tLess(wildTetra,*match)) ? index : -1;
}


// Return index of tetrahedron with all edges shared, to start FindTetra()

G4int G4CMPBiLinearInterp::FirstInteriorTetra() {
  G4int minIndex = Neighbors.size()/4;

  for (G4int i=0; i<(G4int)Neighbors.size(); i++) {
    if (Neighbors[i][0]>minIndex && Neighbors[i][1]>minIndex &&
	Neighbors[i][2]>minIndex) return i;
  }

  return Neighbors.size()/2;
}

// Evaluate mesh at arbitrary location, returning potential or gradient

G4double 
G4CMPBiLinearInterp::GetValue(const G4double pos[2], G4bool quiet) const {
  G4double bary[3] = { 0. };
  FindTetrahedron(&pos[0], bary, quiet);
  staleCache = true;
    
  if (TetraIdx == -1) return 0;
  else {
    return(V[Tetrahedra[TetraIdx][0]] * bary[0] +
           V[Tetrahedra[TetraIdx][1]] * bary[1] +
           V[Tetrahedra[TetraIdx][2]] * bary[2]);    
  }
}

G4ThreeVector 
G4CMPBiLinearInterp::GetGrad(const G4double pos[2], G4bool quiet) const {
  G4double bary[3] = { 0. };
  G4int oldIdx = TetraIdx;
  FindTetrahedron(pos, bary, quiet);

  if (TetraIdx == -1) cachedGrad.set(0.,0.,0.);
  else if (TetraIdx != oldIdx || staleCache) {
    G4double ET[3][2];
    BuildT3x2(ET);
    for (size_t i=0; i<2; ++i) {
      cachedGrad[i] = V[Tetrahedra[TetraIdx][0]]*ET[0][i] +
                      V[Tetrahedra[TetraIdx][1]]*ET[1][i] +
                      V[Tetrahedra[TetraIdx][2]]*ET[2][i];
    }
    cachedGrad[2] = 0.;
    staleCache = false;
  }

  return cachedGrad;
}

void 
G4CMPBiLinearInterp::FindTetrahedron(const G4double pt[2], G4double bary[3],
				      G4bool quiet) const {
  const G4double barySafety = -1e-10;	// Deal with points close to edges

  G4int minBaryIdx = -1;

  G4double bestBary = 0.;	// Norm of barycentric coordinates (below)
  G4int bestTet = -1;

  if (TetraIdx == -1) TetraIdx = TetraStart;

#ifdef G4CMPTLI_DEBUG
  G4cout << "FindTetrahedron pt " << pt[0] << " " << pt[1]
	 << "\n starting from TetraIdx " << TetraIdx << G4endl;
#endif

  for (size_t count = 0; count < Tetrahedra.size(); ++count) {
    Cart2Bary(pt,bary);		// Get barycentric coord in current tetrahedron

#ifdef G4CMPTLI_DEBUG
    G4cout << " Loop " << count << ": Tetra " << TetraIdx << ": "
	   << Tetrahedra[TetraIdx] << "\n bary " << bary[0] << " " << bary[1]
	   << " " << bary[2] << " norm " << BaryNorm(bary)
	   << G4endl;
#endif

    // Point is inside current tetrahedron (TetraIdx)
    if (std::all_of(bary, bary+3,
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
	G4cerr << "G4CMPBiLinearInterp::FindTetrahedron:"
	       << " Point outside of hull!\n pt = "
	       << pt[0] << " " << pt[1] << G4endl;
      }
      return;
    }
  }	// for (size_t count=0 ...

  TetraIdx = bestTet;
  Cart2Bary(pt,bary);
#ifdef G4CMPTLI_DEBUG
  G4cout << "Tetrahedron not found! Using bestTet " << bestTet << " bary "
	 << bary[0] << " " << bary[1] << " " << bary[2]
         << G4endl;
#endif
}

void G4CMPBiLinearInterp::Cart2Bary(const G4double pt[2], G4double bary[3]) const {
  G4double T[2][2];
  G4double invT[2][2];

  for(G4int dim=0; dim<2; ++dim)
    for(G4int vert=0; vert<2; ++vert)
      T[dim][vert] = X[Tetrahedra[TetraIdx][vert]][dim] - X[Tetrahedra[TetraIdx][2]][dim];

  MatInv(T, invT);

  for(G4int k=0; k<2; ++k)
    bary[k] = invT[k][0]*(pt[0] - X[Tetrahedra[TetraIdx][2]][0]) +
              invT[k][1]*(pt[1] - X[Tetrahedra[TetraIdx][2]][1]);

  bary[2] = 1.0 - bary[0] - bary[1] - bary[2];
}

G4double G4CMPBiLinearInterp::BaryNorm(G4double bary[3]) const {
  return (bary[0]*bary[0]+bary[1]*bary[1]+bary[2]*bary[2]);
}

void G4CMPBiLinearInterp::BuildT3x2(G4double ET[3][2]) const {
  G4double T[2][2];
  G4double Tinv[2][2];
  for(G4int dim=0; dim<2; ++dim)
    for(G4int vert=0; vert<2; ++vert)
      T[dim][vert] = (X[Tetrahedra[TetraIdx][vert]][dim]
		      - X[Tetrahedra[TetraIdx][2]][dim]);

  MatInv(T, Tinv);
  for(G4int i=0; i<2; ++i) {
    for(G4int j=0; j<2; ++j)
      ET[i][j] = Tinv[i][j];
    ET[2][i] = -Tinv[0][i] - Tinv[1][i];
  }
}

G4double G4CMPBiLinearInterp::Det2(const G4double matrix[2][2]) const {
  return(matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]);
}

void G4CMPBiLinearInterp::MatInv(const G4double matrix[2][2],
				 G4double result[2][2]) const {
  G4double determ = Det2(matrix);
  if (!(determ == determ)) {
    G4cerr << "WARNING MatInv got NaN determ!"
	   <<"\n " << matrix[0][0] << " " << matrix[0][1]
	   <<"\n " << matrix[1][0] << " " << matrix[1][1]
	   << G4endl;

    for (size_t i=0; i<2; i++) {
      for (size_t j=0; j<2; j++) {
	result[i][j] = 0.;
      }
    }

    return;
  }

  if (fabs(determ) < 1e-9) {
    G4cerr << "WARNING MatInv got determ " << determ << ": zero result from"
	   <<"\n " << matrix[0][0] << " " << matrix[0][1]
	   <<"\n " << matrix[1][0] << " " << matrix[1][1]
	   << G4endl;

    for (size_t i=0; i<2; i++) {
      for (size_t j=0; j<2; j++) {
	result[i][j] = 0.;
      }
    }

    return;
  }

  result[0][0] =  matrix[1][1]/determ;
  result[0][1] = -matrix[0][1]/determ;

  result[1][0] = -matrix[1][0]/determ;
  result[1][1] =  matrix[0][0]/determ;
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
