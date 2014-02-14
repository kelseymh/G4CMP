#include "TriLinearInterp.hh"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <ctime>
#include <map>

using namespace orgQhull;
using std::map;

TriLinearInterp::TriLinearInterp(const vector<vector<G4double> >& xyz, 
				 const vector<G4double>& v)
  : X(xyz), V(v), TetraIdx(0) {
  BuildTetraMesh(xyz);
}

void TriLinearInterp::BuildTetraMesh(const vector<vector<G4double> >& /*xyz*/) {
    time_t start, fin;
    G4cout << "    TriLinearInterp::Constructor: Creating Tetrahedral Mesh..." << G4endl;
    std::time(&start);
    /* Qhull requires a column-major array of the
     * 3D points. i.e., [x1,y1,z1,x2,y2,z2,...]
     */
    double* boxPoints = new double[3*X.size()];
    //G4cout << "Size of input = " << X.size() << " points" << G4endl;
      
    for(size_t i=0; i<X.size(); ++i)
    {
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
     *      d = delaunay triangulation
     *      Qt = triangulate output (makes all regions simplicial)
     *      Qbb = Scales the paraboloid that Qhull creates. This helps with precision.
     *      Qz = Add a point at infinity. This somehow helps with precision...
     */
    Qhull hull = Qhull("", 3, X.size(), boxPoints, "d Qt Qz Qbb");
        
    QhullFacet facet, neighbor;
    QhullVertex vertex;
    QhullLinkedList<QhullFacet>::iterator fItr;
    QhullSet<QhullFacet>::iterator nItr;
    QhullSet<QhullVertex>::iterator vItr;
    map<G4int, G4int> ID2Idx;
    G4int numTet = 0, j;
    vector<vector<G4int> > tmpTetrahedra = vector<vector<G4int> >(hull.facetCount(), vector<G4int>(4, 0));
    vector<vector<G4int> > tmpNeighbors = vector<vector<G4int> >(hull.facetCount(), vector<G4int>(4, -1));
    for (fItr = hull.facetList().begin();fItr != hull.facetList().end(); fItr++)
    {
      facet = *fItr;
      if (!facet.isUpperDelaunay())
      {
        ID2Idx[facet.id()] = numTet;
        j = 0;
        for (vItr = facet.vertices().begin(); vItr != facet.vertices().end(); vItr++)
        {
          vertex = *vItr;
          tmpTetrahedra[numTet][j++] = FindPointID(vertex.point().toStdVector(), vertex.id());
        }
        j = 0;
        for (nItr = facet.neighborFacets().begin(); nItr != facet.neighborFacets().end(); nItr++)
        {
          neighbor = *nItr;
          if (!neighbor.isUpperDelaunay())
            tmpNeighbors[numTet][j++] = neighbor.id();
        }
        ++numTet;
      }
    }
    
//    Tetrahedra.resize(i);
//    Neighbors.resize(i);
    
    for (int i = 0; i < numTet; ++i)
      for (j = 0; j < 4; ++j)
      {
        if (tmpNeighbors[i][j] != -1)
          tmpNeighbors[i][j] = ID2Idx[tmpNeighbors[i][j]];
      }

    tmpTetrahedra.resize(numTet);
    tmpNeighbors.resize(numTet);

    Tetrahedra.swap(tmpTetrahedra);
    Neighbors.swap(tmpNeighbors);

    delete[] boxPoints;		// Clean up local array

    std::time(&fin);
    G4cout << "    TriLinearInterp::Constructor: Took " << difftime(fin, start) << " seconds." << G4endl;
}

G4int TriLinearInterp::FindPointID(const vector<G4double>& point, const G4int id) const
{
  static map<G4int, G4int> qhull2x;
  if (qhull2x.count(id))
  {
    return qhull2x[id];
  }
	
//  for (int i = 0; i < X.size(); ++i)
//  {
//    if (X[i][0] == point[0] && X[i][1] == point[1] && X[i][2] == point[2])
//    {
//      qhull2x[id] = i;
//      return i;
//    }
//  }

  G4int min = 0, max = X.size()-1, mid = X.size()/2;
  G4int d = 0;
  while(true)
  {
    if (X[mid][0] == point[0])
    {
      d = 0;
      while(X[mid+d][0] == point[0] || X[mid-d][0] == point[0])
      {
        if (X[mid+d][0] == point[0] && X[mid+d][1] == point[1] && X[mid+d][2] == point[2])
        {
          qhull2x[id] = mid+d;
          return mid + d;
        }
        else if (X[mid-d][0] == point[0] && X[mid-d][1] == point[1] && X[mid-d][2] == point[2])
        {
          qhull2x[id] = mid-d;
          return mid - d;
        }
        ++d;
      }
      while(X[mid-d][0] == point[0])
      {
        if (X[mid-d][0] == point[0] && X[mid-d][1] == point[1] && X[mid-d][2] == point[2])
        {
          qhull2x[id] = mid-d;
          return mid - d;
        }
        ++d;
      }
    }
    else if (X[mid][0] < point[0])
    {
      min = mid + 1;
      mid = min + ceil((max - min)/2);
    }
    else
    {
      max = mid - 1;
      mid = min + floor((max - min)/2);
    }
  }
}

G4double TriLinearInterp::GetPotential(const G4double pos[3]) const
{
    /* The barycentric coordinates of pos[]*/
    //vector<G4double> bary(4,0);
    G4double bary[4];
    FindTetrahedon(&pos[0], bary);
    
    return(V[Tetrahedra[TetraIdx][0]] * bary[0] +
           V[Tetrahedra[TetraIdx][1]] * bary[1] +
           V[Tetrahedra[TetraIdx][2]] * bary[2] +
           V[Tetrahedra[TetraIdx][3]] * bary[3]);
}

void TriLinearInterp::GetField(const G4double pos[4], G4double field[6]) const
{
    G4double bary[4];
    FindTetrahedon(pos, bary);

    if (TetraIdx == -1)
      for (int i = 0; i < 6; ++i)
        field[i] = 0.0;
    else
    {
      G4double ET[4][3];
      BuildT4x3(pos, ET);
      for (int i = 0; i < 3; ++i)
      {
        field[i] = 0.0;
        field[3+i] = V[Tetrahedra[TetraIdx][0]]*ET[0][i] +
                     V[Tetrahedra[TetraIdx][1]]*ET[1][i] +
                     V[Tetrahedra[TetraIdx][2]]*ET[2][i] +
                     V[Tetrahedra[TetraIdx][3]]*ET[3][i];
      }
    }
}

/*void TriLinearInterp::CalculateTetrahedra()
{
    TetraMesh Triangulation(X);
}*/


void TriLinearInterp::FindTetrahedon(const G4double point[4], G4double bary[4]) const
{
//  G4int minBaryIdx;
  if (TetraIdx == -1) TetraIdx = 0;
  vector<G4int> negBary;
  negBary.reserve(4);
  while (true)
  {
    Cart2Bary(point,bary);

    if (bary[0] > -1.0e-6 && bary[1] > -1.0e-6 && bary[2] > -1.0e-6 && bary[3] > -1.0e-6)
    {
      return;
    }

//    minBaryIdx = 0;
    negBary.clear();
    for (int i = 0; i < 4; ++i)
    {
//      if (bary[i] <= -1.0e-6 && Neighbors[TetraIdx][i] != -1)
      if (bary[i] <= -1.0e-6)
      {
//        minBaryIdx = i;
        negBary.push_back(Neighbors[TetraIdx][i]);
      }
    }
//    TetraIdx = Neighbors[TetraIdx][minBaryIdx];
//    if (TetraIdx == -1)
//    {
//      G4cout << "Point outside of hull!" <<G4endl;
//      TetraIdx = -1;
//      G4cout << point[0]/m << " " << point[1]/m << " " << point[2]/m << G4endl;
//      return;
//    }
      //G4cout << "Before TetraIdx" <<G4endl;
      //G4cout << negBary.size() <<G4endl;
//      G4bool hasInnerNeighbor = false;
//      for (int i = 0; i < negBary.size(); ++i)
//      	hasInnerNeighbor = hasInnerNeighbor || negBary[i] != -1;
      //if (hasInnerNeighbor)

//      if (negBary.size()>0)
        TetraIdx = negBary[round(G4UniformRand()*(negBary.size()-1))];
//      else
//        TetraIdx = -1;
      //G4cout << "After TetraIdx" <<G4endl;
      //G4cout <<TetraIdx<<G4endl;
      if (TetraIdx == -1)
      {
        G4cout << "Outside of hull!" << G4endl;
        G4cout << point[0]/m << " " << point[1]/m << " " << point[2]/m << G4endl;
        return;
      }
  }
}

inline void TriLinearInterp::Cart2Bary(const G4double point[4], G4double bary[4]) const
{
  G4double T[3][3];
  G4double invT[3][3];

  for(int dim=0; dim<3; ++dim)
    for(int vert=0; vert<3; ++vert)
      T[dim][vert] = X[Tetrahedra[TetraIdx][vert]][dim] - X[Tetrahedra[TetraIdx][3]][dim];

  MatInv(T, invT);

  for(int k=0; k<3; ++k)
    bary[k] = invT[k][0]*(point[0] - X[Tetrahedra[TetraIdx][3]][0]) +
              invT[k][1]*(point[1] - X[Tetrahedra[TetraIdx][3]][1]) +
              invT[k][2]*(point[2] - X[Tetrahedra[TetraIdx][3]][2]);

  bary[3] = 1.0 - bary[0] - bary[1] - bary[2];
}

inline void TriLinearInterp::BuildT4x3(const G4double point[4],
				       G4double ET[4][3]) const
{
  G4double T[3][3];
  G4double Tinv[3][3];
  for(int dim=0; dim<3; ++dim)
    for(int vert=0; vert<3; ++vert)
      T[dim][vert] = X[Tetrahedra[TetraIdx][vert]][dim] - X[Tetrahedra[TetraIdx][3]][dim];

  MatInv(T, Tinv);
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
      ET[i][j] = Tinv[i][j];
    ET[3][i] = -Tinv[0][i] - Tinv[1][i] - Tinv[2][i];
  }
}

inline G4double TriLinearInterp::Det3(const G4double matrix[3][3]) const
{
        return(matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])
                -matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2])
                +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]));
}

inline void TriLinearInterp::MatInv(const G4double matrix[3][3], G4double result[3][3]) const
{
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
