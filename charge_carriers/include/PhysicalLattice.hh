#ifndef PhysicalLattice_h
#define PhysicalLattice_h 1

#include "LogicalLattice.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"

class G4VPhysicalVolume;


class PhysicalLattice{

private:
  G4double theta, phi;
  LogicalLattice* Lattice;
  G4VPhysicalVolume* Volume;

  double A;       //Scaling constant for Anh.Dec. mean free path
  double B;       //Scaling constant for Iso.Scat. mean free path
  double dosL;    //Density of states for L-phonons
  double dosST;   //Density of states for ST-phonons
  double dosFT;    //Density of states for FT-phonons
  double beta, gamma, lambda, mu; //dynamical constants for material


public:
  G4AffineTransform LocalToGlobal;
  G4AffineTransform GlobalToLocal;
  
  void setDynamicalConstants(double, double, double, double);
  void setScatteringConstant(G4double);
  void setAnhDecConstant(G4double);
  void setLDOS(double);
  void setSTDOS(double);
  void setFTDOS(double);
  
  double getBeta();
  double getGamma();
  double getLambda();
  double getMu();
  G4double getScatteringConstant();
  G4double getAnhDecConstant();
  double getLDOS();
  double getSTDOS();
  double getFTDOS();


  PhysicalLattice();
  PhysicalLattice(G4VPhysicalVolume*, LogicalLattice*);
  ~PhysicalLattice();
  double mapKtoV(int, G4ThreeVector);
  G4ThreeVector mapKtoVDir(int, G4ThreeVector);
  G4VPhysicalVolume* getVolume();
  void setPhysicalVolume(G4VPhysicalVolume*);
  void setLogicalLattice(LogicalLattice*);
  void setLatticeOrientation(G4double, G4double);
  void setMillerOrientation(int, int, int);


};

#endif
