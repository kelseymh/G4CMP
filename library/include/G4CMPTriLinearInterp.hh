/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$

#ifndef G4CMPTriLinearInterp_h 
#define G4CMPTriLinearInterp_h 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include <vector>
#include <map>
#include <array>

using point = std::array<G4double, 3>;


class G4CMPTriLinearInterp {
public:
  G4CMPTriLinearInterp() : TetraIdx(0) {;}	// Uninitialized version

  G4CMPTriLinearInterp(const std::vector<point >& xyz,
		       const std::vector<G4double>& v);
  ~G4CMPTriLinearInterp() {;}

  // User initialization or re-initialization
  void UseMesh(const std::vector<point >& xyz,
	       const std::vector<G4double>& v);
  
  G4double GetValue(const G4double pos[3]) const;
  G4ThreeVector GetGrad(const G4double pos[3]) const;
  
private:
  std::map<G4int,G4int> qhull2x;
  std::vector<point > X;
  std::vector<G4double> V;
  std::vector<std::array<G4int, 4> > Tetrahedra;
  std::vector<std::array<G4int, 4> > Neighbors;
  mutable G4int TetraIdx;
  mutable G4ThreadLocal G4ThreeVector cachedGrad;
  mutable G4bool staleCache;

  void BuildTetraMesh();	// Builds mesh from pre-initialized 'X' array
  
  void FindTetrahedron(const G4double point[4], G4double bary[4]) const;
  G4int FindPointID(const std::vector<G4double>& point, const G4int id);
  
  void Cart2Bary(const G4double point[4], G4double bary[4]) const;
  void BuildT4x3(G4double ET[4][3]) const;
  G4double Det3(const G4double matrix[3][3]) const;
  void MatInv(const G4double matrix[3][3],
		     G4double result[3][3]) const;
};

#endif
