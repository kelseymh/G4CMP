#include "LogicalLattice.hh"
#include "Randomize.hh"
#include <math.h>
#define PI 3.14159265

LogicalLattice::LogicalLattice(){
  vresTheta=0;
  vresPhi=0;
  dresTheta=0;
  dresPhi=0;
  A=0;
  B=0;
  dosL=0;
  dosST=0;
  dosFT=0;
  beta=0;
  lambda=0;
  gamma=0;
  mu=0;
}

LogicalLattice::~LogicalLattice(){;}


void LogicalLattice::setDynamicalConstants(double Beta,double Gamma, double Lambda, double Mu)
{
  beta=Beta;
  gamma=Gamma;
  lambda=Lambda;
  mu=Mu;

}

void LogicalLattice::setScatteringConstant(G4double b)
{
  B=b;
}

void LogicalLattice::setAnhDecConstant(G4double a)
{
  A=a;
}
void LogicalLattice::setLDOS(double LDOS){
  dosL=LDOS;
}

void LogicalLattice::setSTDOS(double STDOS){
  dosST=STDOS;
}

void LogicalLattice::setFTDOS(double FTDOS){
  dosFT=FTDOS;
}

double LogicalLattice::getBeta()
{
  return beta;
}

double LogicalLattice::getGamma()
{
  return gamma;
}

double LogicalLattice::getLambda()
{
  return lambda;
}

double LogicalLattice::getMu()
{
  return mu;
}

G4double LogicalLattice::getScatteringConstant()
{
  return B;
}

G4double LogicalLattice::getAnhDecConstant()
{
  return A;
}

double LogicalLattice::getLDOS()
{
  return dosL;
}

double LogicalLattice::getSTDOS()
{
  return dosST;
}

double LogicalLattice::getFTDOS()
{
  return dosFT;
}

///////////////////////////////////////////
//Load map of group velocity scalars (m/s)
////////////////////////////////////////////
bool LogicalLattice::loadMap(int tRes, int pRes, int polarizationState, string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V map exceeds maximum resolution of "<<MAXRES<<" by "<<MAXRES<<". terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  mapFile.clear();

  thetaRes=tRes;
  phiRes=pRes;
  mapFile.open(m.data());
  if(!mapFile.is_open()) return false;
  
  for(int theta = 0; theta<thetaRes; theta++){
    for(int phi = 0; phi<phiRes; phi++){
	mapFile>>map[polarizationState][theta][phi];	
	//G4cout<<"\nLogicalLattice::loadMap:loading map "<<m<<" theta: "<<theta<<" phi: "<<phi<<" map: "<<map[polarizationState][theta][phi];
    }
  }
  mapFile.close();
  G4cout<<"\nLogicalLattice::loadMap() sucessful (Vg scalars)."<<endl;
  vresTheta=tRes; //store map dimensions
  vresPhi=pRes;
  return true;
}

////////////////////////////////////
//Load map of group velocity unit vectors
///////////////////////////////////
bool LogicalLattice::load_NMap(int tRes, int pRes, int polarizationState, string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V map exceeds maximum resolution of "<<MAXRES<<" by "<<MAXRES<<" terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  mapFile.clear();

  thetaRes=tRes;
  phiRes=pRes;
  mapFile.open(m.data());
  if(!mapFile.is_open()) return false;
  
  for(int theta = 0; theta<thetaRes; theta++){
    for(int phi = 0; phi<phiRes; phi++){
      for(int coord = 0; coord<3; coord++){
	mapFile>>n_map[polarizationState][theta][phi][coord];
      }
      //G4cout<<"\n theta: " <<theta<<" phi: "<<phi<<" [xyz]: "<<n_map[polarizationState][theta][phi][0]<<" "<<n_map[polarizationState][theta][phi][1]<<" "<<n_map[polarizationState][theta][phi][2];
    }
  }
  G4cout<<"\n";
  mapFile.close();
  G4cout<<"\nLogicalLattice::load_NMap() sucessful"<<endl;
  dresTheta=tRes; //store map dimesnions
  dresPhi=pRes;
  return true;
}


double LogicalLattice::mapKtoV(int polarizationState,G4ThreeVector k)
{
  double theta, phi, tRes, pRes;

  tRes=PI/(vresTheta);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*PI/(vresPhi);

  
  theta=k.getTheta();
  phi=k.getPhi();

  if(phi<0) phi = phi + 2*PI;
  if(theta>PI) theta=theta-PI;
  //phi=[0 to 2 PI] in accordance with DMC //if(phi>PI/2) phi=phi-PI/2;
  if(map[polarizationState][int(theta/tRes)][int(phi/pRes)]==0){
      G4cout<<"\nFound v=0 for polarization "<<polarizationState<<" theta "<<theta<<" phi "<<phi<< " translating to map coords " << "theta "<< int(theta/tRes) << " phi " << int(phi/pRes)<<endl;
  }

  //  G4cout<<"\nTheta: "<<theta<<" theta index: "<<int(theta/tRes)<<" tRes: "<<tRes;
  return map[polarizationState][int(theta/tRes)][int(phi/pRes)];
  
  ///////////////////debugging purposes///////////////////////////
  //return 5000.0;
}


G4ThreeVector LogicalLattice::mapKtoVDir(int polarizationState,G4ThreeVector k)
{  
  double theta, phi, tRes, pRes;

  tRes=PI/(dresTheta-1);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*PI/(dresPhi-1);

  theta=k.getTheta();
  phi=k.getPhi(); 

  if(theta>PI) theta=theta-PI;
  //phi=[0 to 2 PI] in accordance with DMC //if(phi>PI/2) phi=phi-PI/2;
  if(phi<0) phi = phi + 2*PI;

  G4int iTheta = int(theta/tRes+0.5);
  G4int iPhi = int(phi/pRes+0.5);

  
  G4ThreeVector v(n_map[polarizationState][iTheta][iPhi][0],
		  n_map[polarizationState][iTheta][iPhi][1],
		  n_map[polarizationState][iTheta][iPhi][2]);


  //////debugging purposes only//////////////
  //v.set(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  //v.set(k.getX(), k.getY(), k.getZ());

  return v.unit();
}
