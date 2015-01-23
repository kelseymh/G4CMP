// $Id$

#ifndef G4CMPTriLinearInterp_h 
#define G4CMPTriLinearInterp_h 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>


class G4CMPTriLinearInterp {
public:
  G4CMPTriLinearInterp() : TetraIdx(0) {;}	// Uninitialized version

  G4CMPTriLinearInterp(const std::vector<std::vector<G4double> >& xyz,
		       const std::vector<G4double>& v);
  ~G4CMPTriLinearInterp() {;}

  // User initialization or re-initialization
  void UseMesh(const std::vector<std::vector<G4double> >& xyz,
	       const std::vector<G4double>& v);
  
  G4double GetPotential(const G4double pos[3]) const;
  void GetField(const G4double pos[4], G4double field[6]) const;
  
private:
  std::vector<std::vector<G4double> > X;
  std::vector<G4double> V;
  std::vector<std::vector<G4int> > Tetrahedra;
  std::vector<std::vector<G4int> > Neighbors;
  mutable G4int TetraIdx;

  void BuildTetraMesh();	// Builds mesh from pre-initialized 'X' array
  
  void FindTetrahedron(const G4double point[4], G4double bary[4]) const;
  G4int FindPointID(const std::vector<G4double>& point, const G4int id) const;
  
  void Cart2Bary(const G4double point[4], G4double bary[4]) const;
  void BuildT4x3(G4double ET[4][3]) const;
  G4double Det3(const G4double matrix[3][3]) const;
  void MatInv(const G4double matrix[3][3],
		     G4double result[3][3]) const;
};

#endif
