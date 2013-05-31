#ifndef LogicalLattice_h
#define LogicalLattice_h

#define MAXRES 322  //maximum one dimensional map resolution

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"


using namespace std;

class LogicalLattice{

private:
  
  int vresTheta; //velocity  map theta resolution (inclination)
  int vresPhi;   //velocity  map phi resolution  (azimuth)
  int dresTheta; //direction map theta resn
  int dresPhi;   //direction map phi resn 

  double map[3][MAXRES][MAXRES];  //map for group velocity scalars
  double n_map[3][MAXRES][MAXRES][3];//normalized map containing group V direction unit vectors

  double A;       //Scaling constant for Anh.Dec. mean free path
  double B;       //Scaling constant for Iso.Scat. mean free path
  double dosL;    //Density of states for L-phonons
  double dosST;   //Density of states for ST-phonons
  double dosFT;    //Density of states for FT-phonons
  double beta, gamma, lambda, mu; //dynamical constants for material


  ifstream mapFile;
  int thetaRes, phiRes;


public:


  
  LogicalLattice();
  ~LogicalLattice();

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
  
  bool loadMap(int, int, int, string);
  bool load_NMap(int, int, int, string);
  double mapKtoV(int, G4ThreeVector);   //get full group velocity vector
  G4ThreeVector mapKtoVDir(int, G4ThreeVector);//get normalized group velocity direction so that normalisatioon does not have to be done at run time

};

#endif
