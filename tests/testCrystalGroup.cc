#include "G4CMPCrystalGroup.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <cmath>

using std::cout; 
using std::endl;

G4double xyangle(const G4CMPCrystalGroup& crystal) {
  return std::acos(crystal.axis[0].dot(crystal.axis[1]));
}

G4double xzangle(const G4CMPCrystalGroup& crystal) {
  return std::acos(crystal.axis[0].dot(crystal.axis[2]));
}

G4double yzangle(const G4CMPCrystalGroup& crystal) {
  return std::acos(crystal.axis[1].dot(crystal.axis[2]));
}

int main() {
  cout << "Testing G4CMPCrystalGroup for non-orthogonal axes" << endl;

  G4CMPCrystalGroup rhomb(G4CMPCrystalGroup::rhombohedral, 60.*deg);

  G4cout << rhomb.Name() << " created w/60 deg angles :"
	 << " X-Y " << xyangle(rhomb)/deg << " X-Z " << xzangle(rhomb)/deg
	 << " Y-Z " << yzangle(rhomb)/deg << G4endl;

  G4CMPCrystalGroup tricl(60.*deg, 70.*deg, 80.*deg);	// Triclinic

  G4cout << tricl.Name() << " created with 60, 70, 80 deg angles :"
	 << " X-Y " << xyangle(tricl)/deg << " X-Z " << xzangle(tricl)/deg
	 << " Y-Z " << yzangle(tricl)/deg << G4endl;

  rhomb.Set(G4CMPCrystalGroup::cubic);

  G4cout << rhomb.Name() << " has angles :"
	 << " X-Y " << xyangle(rhomb)/deg << " X-Z " << xzangle(rhomb)/deg
	 << " Y-Z " << yzangle(rhomb)/deg << G4endl;

}
