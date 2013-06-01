#ifndef Interpolation3D_h
#define Interpolation3D_h
#include"globals.hh"
#include<vector>

using std::vector;
class Interpolation3D
{
    public:
        Interpolation3D();
        Interpolation3D(const vector<vector<G4double> >& xyz, 
                        const vector<G4double>& v);
        ~Interpolation3D() {}
        G4double GetValue(const G4double pos[3]);
    private:
        vector<vector<G4double> > X;
        vector<G4double> V;
        vector<vector<G4int> > tetrahedra;
        void CalculateTetrahedra();
        G4double Det(const vector<vector<G4double> >& matrix);
        void MatInv(const vector<vector<G4double> >& matrix, 
                vector<vector<G4double> >& result);
        void Cart2Bary(const vector<vector<G4double> >& xyzMat, 
                const vector<vector<G4int> >& tetras, const G4double point[3], 
                G4int& tetraIdx, vector<G4double>& bary);
};
#endif
