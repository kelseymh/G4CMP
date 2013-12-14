#ifndef Triangulation3D_h
#define Triangulation3D_h
#include<globals.hh>
#include<vector>

using std::vector;
class Triangulation3D
{
    public:
        Triangulation3D(const vector<vector<G4double> >& xyz);
        ~Triangulation3D();
        inline void GetTetrahedra(vector<vector<G4int> >& tetra){tetra=tetraIndices;}
        inline void GetXYZMatrix(vector<vector<G4double> >& mat){mat=X;}
    private:
        vector<vector<G4double> > X;
        vector<vector<G4int> > tetraIndices;
        void RunQhull();
        G4double Det(const vector<vector<G4double> >& matrix);
};
#endif
