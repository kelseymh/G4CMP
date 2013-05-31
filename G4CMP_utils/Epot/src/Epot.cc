#include <fstream>
#include <iostream>
#include <math.h>
#include "float.h" //For DBL_EPSILON, DBL_MAX, etc
#include "Epot.hh"

using std::string;
using std::vector;
using std::cout;
using std::endl;

DetectorField::DetectorField(G4String filename): xmin(0),xmax(0),
                                                 ymin(0),ymax(0),
                                                 zmin(0),zmax(0)
{
    // Read in File that has x,y,z,V
    std::fstream input(filename.c_str());
	    
    G4String buffer;
    for(int i=0; i<8; ++i)
        getline(input, buffer);

    double x,y,z,v;
    int lineNum=0;
    int n=0;
    vector<double> tempVec(3,0);
    while(!input.eof())
    {
        input >> x >> y >> z >> v;

        if(x < xmin)
            xmin = x;
        else if(x > xmax)
            xmax = x;
        if(y < ymin)
            ymin = y;
        else if(y > ymax)
            ymax = y;
        if(z < zmin)
            zmin = z;
        else if(z > zmax)
            zmax = z;

        if(xyzMat.size()>=100000*n)
        {
            ++n;
            xyzMat.reserve(100000*n);
            V.reserve(100000*n);
        }
            
        tempVec[0] = x;
        tempVec[1] = y;
        tempVec[2] = z;
        xyzMat.push_back(tempVec);
        V.push_back(v);

        ++lineNum;
    }
    xyzMat.resize(lineNum);
    V.resize(lineNum);

    input.close();
    Interpolation = Interpolation3D(xyzMat, V);
}

DetectorField::~DetectorField()
{
}

void DetectorField::GetFieldValue(const G4double Point[4], G4double *Efield) const
{
    /*
     * To calculate the field value, we must first calculate the
     * voltage at several points and then calulate the gradient.
     * To calculate the voltages we will use the Qhull tetrahedral mesh.
     */

    /* Strip off the units to be consistent with Qhull */
    double newPoint[4] = {Point[0]/m, Point[1]/m, Point[2]/m, Point[3]/m};

    if((newPoint[0]<xmin || newPoint[1]>xmax) ||
       (newPoint[1]<ymin || newPoint[1]>ymax) ||
       (newPoint[2]<zmin || newPoint[2]>zmax))
        G4cout << "Warning: The point " << newPoint[0] << "m, "
               << newPoint[1] << "m , " << newPoint[2]
               << "m is clearly outside of the cubic volume defined by\n"
               << "         the minimum and maximum values of x,y,z from"
               << " the input file."
               << G4endl;
        
    /* https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic */
    G4double dx = sqrt(DBL_EPSILON)*newPoint[0];
    G4double dy = sqrt(DBL_EPSILON)*newPoint[1];
    G4double dz = sqrt(DBL_EPSILON)*newPoint[2];

    if(dx == 0)
        dx = sqrt(DBL_EPSILON);
    if(dy == 0)
        dy = sqrt(DBL_EPSILON);
    if(dz == 0)
        dz = sqrt(DBL_EPSILON);

    /* Create x1 and x2 for deltaX in deltaV/deltaX */
    G4double x1 = newPoint[0] - dx;
    G4double x2 = newPoint[0] + dx;
    G4double y1 = newPoint[1] - dy;
    G4double y2 = newPoint[1] + dy;  
    G4double z1 = newPoint[2] - dz;
    G4double z2 = newPoint[2] + dz;

    /* deltaV = V(gradPointX2) - V(gradPointX1) */
    G4double gradPointX1[3] = {x1,newPoint[1],newPoint[2]};
    G4double gradPointX2[3] = {x2,newPoint[1],newPoint[2]};
    G4double gradPointY1[3] = {newPoint[0],y1,newPoint[2]};
    G4double gradPointY2[3] = {newPoint[0],y2,newPoint[2]};
    G4double gradPointZ1[3] = {newPoint[0],newPoint[1],z1};
    G4double gradPointZ2[3] = {newPoint[0],newPoint[1],z2};

    G4double V_X1, V_X2, V_Y1, V_Y2, V_Z1, V_Z2;

    Interpolation3D pnt = Interpolation;
    V_X1 = pnt.GetValue(gradPointX1);
    V_X2 = pnt.GetValue(gradPointX2);
    V_Y1 = pnt.GetValue(gradPointY1);
    V_Y2 = pnt.GetValue(gradPointY2);
    V_Z1 = pnt.GetValue(gradPointZ1);
    V_Z2 = pnt.GetValue(gradPointZ2);

    /* Put units back on */
    Efield[0] = -1*(V_X2-V_X1)/(2*dx)*volt/m;
    Efield[1] = -1*(V_Y2-V_Y1)/(2*dy)*volt/m;
    Efield[2] = -1*(V_Z2-V_Z1)/(2*dz)*volt/m;
}
