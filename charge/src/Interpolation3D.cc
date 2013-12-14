#include "Interpolation3D.hh"
#include "Triangulation3D.hh"
#include "libqhullcpp/Qhull.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

#include <ctime>

Interpolation3D::Interpolation3D()
{
    /* Don't use the default constructor. This is here so that classes can
     * have an interpolation object.
     */
    vector<vector<G4double> > tempX(1,vector<G4double>(3,0));
    X.swap(tempX);
    vector<G4double> tempV(1,0);
    V.swap(tempV);
    vector<vector<G4int> > tempTetra(1,vector<G4int>(3,0));
    tetrahedra.swap(tempTetra);
}
Interpolation3D::Interpolation3D(const vector<vector<G4double> >& xyz, 
                                 const vector<G4double>& v) : X(xyz), V(v)
{
    CalculateTetrahedra();
}

G4double Interpolation3D::GetValue(const G4double pos[3])
{
    /* Indices indicating which tetrahedron pos[] is in */
    G4int tetraIdx;
        
    /* The barycentric coordinates of pos[]*/
    vector<G4double> bary(4,0);
    Cart2Bary(X, tetrahedra, &pos[0], tetraIdx, bary);
    
    return(V[tetrahedra[tetraIdx][0]] * bary[0] +
           V[tetrahedra[tetraIdx][1]] * bary[1] +
           V[tetrahedra[tetraIdx][2]] * bary[2] +
           V[tetrahedra[tetraIdx][3]] * bary[3]);
}

void Interpolation3D::CalculateTetrahedra()
{
    Triangulation3D tri(X);
    tri.GetTetrahedra(tetrahedra);
}

void Interpolation3D::Cart2Bary(const vector<vector<G4double> >& xyzMat, const vector<vector<G4int> >& tetras, const G4double point[3], G4int& tetraIdx, vector<G4double>& bary)
{
    /* Search through every tetrahedron until we find the one that point[] is in.
     * We know that is the case when all barycentric coordinates are between
     * 0 and 1. The way they are calculated here normalizes them to 1, so we only
     * have to check if bary[0-3] > 0
     *
     * If the point doesn't appear to be inside a tetrahedron, we find the nearest
     * tetrahedron and incrementally move the point until it's inside that tetrahedron.
     * This does not seem to happen very often at all.
     *
     * This function is very messy and ugly. I'd like to clean it up, both for
     * performance and readability.
     */
 
    //Distance from the point to the current simplex. Used to find nearest simplex. 
    double dist = DBL_MAX; 
    //Did we find the simplex containing the point?
    bool found = false; 

    //Temporary variables are used inside the for-loop.
    vector<double> tmpBary(4,0);

    /*T is the transformation matrix for going between barycentric and cartesian.
    * invT is its inverse.
    */
    vector<vector<double> > T(3,vector<double>(3,0));
    vector<vector<double> > invT(3,vector<double>(3,0));
    int max = tetras.size(); //For the for-loop iterator.
    //time_t start = time(NULL);
    //time_t end;
    for(int tmpTetraIdx=0; tmpTetraIdx<max; ++tmpTetraIdx)
    {
        /* Need to build T matrix, then invert it to calculate barycentric
         * coordinates. For more info, 
         * https://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)#Barycentric_coordinates_on_tetrahedra
         */
        for(int dim=0; dim<3; ++dim)
            for(int vert=0; vert<3; ++vert)
                T[dim][vert] = xyzMat[tetras[tmpTetraIdx][vert]][dim] - 
                               xyzMat[tetras[tmpTetraIdx][3]][dim];

        MatInv(T, invT);

        for(int k=0; k<3; ++k)
            tmpBary[k] = invT[k][0]*(point[0] - xyzMat[tetras[tmpTetraIdx][3]][0]) +
                         invT[k][1]*(point[1] - xyzMat[tetras[tmpTetraIdx][3]][1]) +
                         invT[k][2]*(point[2] - xyzMat[tetras[tmpTetraIdx][3]][2]);

        tmpBary[3] = 1.0 - tmpBary[0] - tmpBary[1] - tmpBary[2];

        if((tmpBary[0] >= 0) && (tmpBary[1] >= 0) && (tmpBary[2] >= 0) &&
           (tmpBary[3] >= 0) && (tmpBary[0] <= 1) && (tmpBary[1] <= 1) &&
           (tmpBary[2] <= 1) && (tmpBary[3] <= 1))
        {
            tetraIdx = tmpTetraIdx;
            bary = tmpBary;
            found = true;
            break;
        }
        else
        {
            if(dist > pow(tmpBary[0] - .25, 2) + pow(tmpBary[1] - .25, 2) + 
                      pow(tmpBary[2] - .25, 2) + pow(tmpBary[3] - .25, 2))
            {
                dist = pow(tmpBary[0] - .25, 2) + pow(tmpBary[1] - .25, 2) + 
                       pow(tmpBary[2] - .25, 2) + pow(tmpBary[3] - .25, 2);
                tetraIdx = tmpTetraIdx;
                bary = tmpBary;
            }
        }
    }
    if(!found)
    {
        //time_t start = time(NULL);
        //time_t end;
        G4cout << "Warning: This point could not be placed in a tetrahedron. \n"
               << "         It may be outside the volume or in a small region that is\n"
               << "         not covered by the tetrahedral mesh. We will joggle the point\n"
               << "         to the nearest tetrahedron and use that value." << G4endl;
        /* Since this point isn't in a tetrahedron, we will iteratively move it
         * toward the center of the nearest tetrahedron. Eventually, it will be
         * inside that tetrahedron.
         */
        std::ofstream outfile("badpoints", std::ios::app);
        outfile << point[0] << " " << point[1] << " " << point[2] << "\n";
        outfile.close();

        //NOTE: tmpBary == bary

        for(int i=0; i<1000; ++i) //Move point toward nearest tetrahedron
        {
            for(int j=0; j<4; ++j)
                tmpBary[j] += (.25 - bary[j])/1000.0;

            if((tmpBary[0] >= 0) && (tmpBary[1] >= 0) && (tmpBary[2] >= 0) &&
               (tmpBary[3] >= 0) && (tmpBary[0] <= 1) && (tmpBary[1] <= 1) &&
               (tmpBary[2] <= 1) && (tmpBary[3] <= 1))
            {
                bary = tmpBary;
                break;
            }
        }
    }
    /*********** Qhull has its own method to find the best tetra, but I don't know if
     *********** it's better yet... So just ignore this stuff
    coordT qh_point[4]={0,0,0,0};
    boolT isoutside;
    realT bestdist;
    facetT *facet1;
    vertexT *vertex, **vertexp;

    for(int i=0; i<3; i++)
    {
        qh_point[i] = point[i];
        qh_point[3] += qh_point[i]*qh_point[i];
    }

    G4cout << "before setdealunay" << G4endl;
    //qh_setdelaunay (4, 1, qh_point);
    G4cout << "After setdelaunay, before bestfacet." << G4endl;
    facet1= qh_findbestfacet (qh_point, qh_ALL, &bestdist, &isoutside);
    vector<vector<double> > T(3,vector<double>(3,0));
    vector<vector<double> > invT(3,vector<double>(3,0));
    int point_id[4];

    int i=0;
    FOREACHvertex_ (facet1->vertices)
    {
        point_id[i++] = qh_pointid(vertex->point);
    }

    for(int dim=0; dim<3; dim++)
        for(int vert=0; vert<3; vert++)
            T[dim][vert] = xyzMat[point_id[vert]][dim] - xyzMat[point_id[3]][dim];

    MatInv(T,invT);

    for(int k=0; k<3; k++)
        bary[k] = invT[k][0]*(point[0] - xyzMat[point_id[3]][0])+
                  invT[k][1]*(point[1] - xyzMat[point_id[3]][1])+
                  invT[k][2]*(point[2] - xyzMat[point_id[3]][2]);

    bary[3] = 1.0 - bary[0] - bary[1] - bary[2];
    */
}

G4double Interpolation3D::Det(const vector<vector<G4double> >& matrix)
{
    int rows = matrix.size(); //rows = columns

    if(rows==2)
        return(matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1]);
    else if(rows==3)
        return(matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])
                -matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2])
                +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]));
    else
    {    
        /*This function is never called for anything but 3x3 and 2x2, but just to 
         * future-proof it...
         */
        vector<vector<double> > newMatrix(rows-1,vector<double>(rows-1,0));

        double determ = 0;
        for(int i=0; i<rows; ++i) //columns
        {
            for(int j=0; j<rows-1; ++j)
                for(int k=0; k<rows-1; ++k) //columns
                    if(k<i)
                        newMatrix[j][k] = matrix[j+1][k];
                    else
                        newMatrix[j][k] = matrix[j+1][k+1];
            if((i+1)%2==0)
                determ += -1*matrix[0][i]*Det(newMatrix);
            else
                determ += matrix[0][i]*Det(newMatrix);
        }
        return(determ);
    }
}

void Interpolation3D::MatInv(const vector<vector<G4double> >& matrix, vector<vector<G4double> >& result)
{
    G4double determ = Det(matrix);
    result[0][0] = (matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])/determ;
    result[1][0] = (matrix[1][2]*matrix[2][0] - matrix[1][0]*matrix[2][2])/determ;
    result[2][0] = (matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0])/determ;

    result[0][1] = (matrix[0][2]*matrix[2][1] - matrix[0][1]*matrix[2][2])/determ;
    result[1][1] = (matrix[0][0]*matrix[2][2] - matrix[2][0]*matrix[0][2])/determ;
    result[2][1] = (matrix[0][1]*matrix[2][0] - matrix[0][0]*matrix[2][1])/determ;

    result[0][2] = (matrix[0][1]*matrix[1][2] - matrix[1][1]*matrix[0][2])/determ;
    result[1][2] = (matrix[1][0]*matrix[0][2] - matrix[0][0]*matrix[1][2])/determ;
    result[2][2] = (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1])/determ;
    /*
    vector<vector<double> > minorMat(2,vector<double>(2,0));

    int row = 0;
    int column = 0;

    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        {
            row = 0;
            for(int l=0; l<3; ++l)
            {
                column = 0;
                if(l != i)
                {
                    for(int k=0; k<3; ++k)
                        if(k != j)
                            minorMat[row][column++] = matrix[k][l]; //matrix is transposed
                    ++row;
                }
            }
            if((i+j)%2==1)
                result[i][j] = -1*1.0/Det(matrix)*Det(minorMat);
            else
                result[i][j] = 1.0/Det(matrix)*Det(minorMat);
        }
        */
}
