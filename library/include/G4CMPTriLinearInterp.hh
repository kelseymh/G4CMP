#ifndef G4CMPTriLinearInterp_h 
#define G4CMPTriLinearInterp_h 

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
//#include "Tetrahedron.hh"

using std::vector;
class G4CMPTriLinearInterp {
  public:
    G4CMPTriLinearInterp( const vector<vector<G4double> >& xyz,
                         const vector<G4double>& v );
    ~G4CMPTriLinearInterp() {}

    G4double GetPotential( const G4double pos[3] ) const;
    void GetField( const G4double pos[4], G4double field[6] ) const;
        
  private:
    vector<vector<G4double> > X;
    vector<G4double> V;
    vector<vector<G4int> > Tetrahedra;
    vector<vector<G4int> > Neighbors;
    mutable G4int TetraIdx;

    void BuildTetraMesh( const vector<vector<G4double> >& xyz );
    void FindTetrahedon( const G4double point[4], G4double bary[4] ) const;
    G4int FindPointID( const vector<G4double>& point, const G4int id ) const;

    inline void Cart2Bary( const G4double point[4], G4double bary[4] ) const;
    inline void BuildT4x3( const G4double point[4], G4double T[4][3] ) const;
    inline G4double Det3( const G4double matrix[3][3] ) const;
    inline void MatInv( const G4double matrix[3][3],
			G4double result[3][3] ) const;
};

#endif
