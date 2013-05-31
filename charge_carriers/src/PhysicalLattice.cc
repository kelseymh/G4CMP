#include "PhysicalLattice.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>

#define PI 3.14159265

PhysicalLattice::PhysicalLattice(){
  Lattice=NULL;
  Volume=NULL;
  theta=0;
  phi=0;
}

PhysicalLattice::PhysicalLattice(G4VPhysicalVolume* Vol, LogicalLattice* Lat){
  Lattice=Lat;
  Volume=Vol;
  A=Lattice->getAnhDecConstant();
  B=Lattice->getScatteringConstant();
  dosL=Lattice->getLDOS();
  dosST=Lattice->getSTDOS();
  dosFT=Lattice->getFTDOS();
  beta=Lattice->getBeta();
  gamma=Lattice->getGamma();
  lambda=Lattice->getLambda();
  mu=Lattice->getMu();

  G4RotationMatrix *rot = Volume->GetObjectRotation();

  GlobalToLocal = G4AffineTransform(*rot);
  LocalToGlobal = GlobalToLocal.Invert();
}

PhysicalLattice::~PhysicalLattice(){;}

void PhysicalLattice::setDynamicalConstants(double Beta, double Gamma, double Lambda, double Mu)
{
  beta=Beta;
  gamma=Gamma;
  lambda=Lambda;
  mu=Mu;
}

void PhysicalLattice::setScatteringConstant(G4double a){
  A=a;
}

void PhysicalLattice::setAnhDecConstant(G4double b){
  B=b;
}

void PhysicalLattice::setLDOS(double LDOS){
  dosL=LDOS;
}
void PhysicalLattice::setSTDOS(double STDOS)
{
  dosST = STDOS;
}
void PhysicalLattice::setFTDOS(double FTDOS){
  dosFT = FTDOS;
}

double PhysicalLattice::getBeta(){
  return beta;
}

double PhysicalLattice::getGamma(){
  return gamma;
}
double PhysicalLattice::getLambda(){
  return lambda;
}

double PhysicalLattice::getMu()
{
  return mu;
}

G4double PhysicalLattice::getScatteringConstant()
{
  return B;
}

G4double PhysicalLattice::getAnhDecConstant()
{
  return A;
}

double PhysicalLattice::getLDOS(){
  return dosL;
}

double PhysicalLattice::getSTDOS(){
  return dosST;
}

double PhysicalLattice::getFTDOS()
{
  return dosFT;
}


///////////////////////////////
//Loads the group velocity in m/s
/////////////////////////////
double PhysicalLattice::mapKtoV(int polarizationState, G4ThreeVector k){
  double groupVelocity;

  k.rotate(G4ThreeVector(0,1,0), theta).rotate(G4ThreeVector(0,0,1), phi);
  groupVelocity = Lattice->mapKtoV(polarizationState, k);
  k.rotate(G4ThreeVector(0,0,1), -phi).rotate(G4ThreeVector(0,1,0), -theta);

  return groupVelocity;
}


///////////////////////////////
//Loads the normalized direction vector along VG
///////////////////////////////
G4ThreeVector PhysicalLattice::mapKtoVDir(int polarizationState, G4ThreeVector k){

  G4ThreeVector GroupVelocity;

  k=k.rotate(G4ThreeVector(0,1,0), theta).rotate(G4ThreeVector(0,0,1), phi);
  GroupVelocity = Lattice->mapKtoVDir(polarizationState, k);
  //return GroupVelocity;
  
  return GroupVelocity.rotate(G4ThreeVector(0,0,1), -phi).rotate(G4ThreeVector(0,1,0), -theta).unit();
  //  return G4ThreeVector(1,0,0).unit();
  //return G4ThreeVector(0,1,0).unit();
}

G4VPhysicalVolume* PhysicalLattice::getVolume(){
  return Volume;
}

void PhysicalLattice::setPhysicalVolume(G4VPhysicalVolume* Vol){
  Volume=Vol;
}

void PhysicalLattice::setLatticeOrientation(G4double t_rot, G4double p_rot){
  theta=t_rot;
  phi= p_rot;
}

void PhysicalLattice::setMillerOrientation(int l, int m, int n){
  theta=PI/2-atan2(n+0.000001,l+0.000001)*rad;
  phi=PI/2-atan2(l+0.000001,m+0.000001)*rad;
}

void PhysicalLattice::setLogicalLattice(LogicalLattice* Lat){
  Lattice=Lat;
}
