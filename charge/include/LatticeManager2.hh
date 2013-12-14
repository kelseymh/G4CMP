#ifndef LatticeManager2_h
#define LatticeManager2_h 1


#define MAXLAT 10 //maximum number of different crystal lattices supported

#include "G4ThreeVector.hh"
#include "PhysicalLattice.hh"

class LatticeManager2
{
private:

 

public:
  LatticeManager2();
  ~LatticeManager2();

  static PhysicalLattice* LatticeList[MAXLAT];
  static int totalLattices;


  static PhysicalLattice* getPhysicalLattice(G4VPhysicalVolume*);
  static bool registerLattice(PhysicalLattice*);
  static bool hasLattice(G4VPhysicalVolume*);
  static double mapKtoV(G4VPhysicalVolume*,int,const G4ThreeVector &);
  static G4ThreeVector mapKtoVDir(G4VPhysicalVolume*,int,const G4ThreeVector &);

};



#endif
