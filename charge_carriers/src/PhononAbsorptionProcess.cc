//#define BETA -0.732
//#define GAMMA -0.708
//#define LAMBDA 0.376
//#define MU 0.561
#define PI 3.1417

#include "PhononAbsorptionProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "Phonon.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "LPhonon.hh"
#include <math.h>

#include "PhononTrackInformation.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "LatticeManager2.hh"



PhononAbsorptionProcess::PhononAbsorptionProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
   Lattice=0;

}

PhononAbsorptionProcess::~PhononAbsorptionProcess()
{

}

PhononAbsorptionProcess::PhononAbsorptionProcess(PhononAbsorptionProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  PhononAbsorptionProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{

  Lattice = LatticeManager2::getPhysicalLattice(aTrack.GetVolume());

  G4double A=Lattice->getAnhDecConstant();
  G4double h=6.626068e-34*m2*kg/s;
  G4double E= aTrack.GetKineticEnergy();
  
  //Calculate mean free path for anh. decay
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*(E/h)*A)*aTrack.GetVelocity();
  //G4cout<<"\nPhononAbsorbtionProces::GetMeanFreePath:Absorption process: Phonon of energy E="<<E/eV<<" and velocity v="<<aTrack.GetVelocity()<<" has mean free path l="<<mfp<<endl;
  
  
  *condition = NotForced;
  return mfp;
}



G4VParticleChange*
  PhononAbsorptionProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&)
{
  
  aParticleChange.Initialize(aTrack);

  Lattice = LatticeManager2::getPhysicalLattice(aTrack.GetVolume());
  BETA=Lattice->getBeta();
  GAMMA=Lattice->getGamma();
  LAMBDA=Lattice->getLambda();
  MU=Lattice->getMu();

  //Destroy the parent phonon and create the daughter phonons.
  //74% chance that daughter phonons are both transverse, 
  //26% Transverse and Longitudinal
  if(G4UniformRand()>0.740) makeLTSecondaries(aTrack); else makeTTSecondaries(aTrack);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);    
       
  return &aParticleChange;
}


G4bool PhononAbsorptionProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  //Only L-phonons decay
  return (&aPD==LPhonon::PhononDefinition());
}

inline double PhononAbsorptionProcess::getLTDecayProb(double d, double x){
  //probability density of energy distribution of L'-phonon in L->L'+T process
  //d=delta= ratio of group velocities vl/vt and x is the fraction of energy in the longitudinal mode, i.e. x=EL'/EL
  return (1/(x*x))*(1-x*x)*(1-x*x)*((1+x)*(1+x)-d*d*((1-x)*(1-x)))*(1+x*x-d*d*(1-x)*(1-x))*(1+x*x-d*d*(1-x)*(1-x));
}

inline double PhononAbsorptionProcess::getTTDecayProb(double d, double x){  
  //probability density of energy distribution of T-phonon in L->T+T process
  
  //dynamic constants from Tamura, PRL31, 1985
  G4double A = 0.5*(1-d*d)*(BETA+LAMBDA+(1+d*d)*(GAMMA+MU));
  G4double B = BETA+LAMBDA+2*d*d*(GAMMA+MU);
  G4double C = BETA + LAMBDA + 2*(GAMMA+MU);
  G4double D = (1-d*d)*(2*BETA+4*GAMMA+LAMBDA+3*MU);

  return (A+B*d*x-B*x*x)*(A+B*d*x-B*x*x)+(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)))*(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)));
}

inline double PhononAbsorptionProcess::makeLDeviation(double d, double x){
  //change in L'-phonon propagation direction after decay

  return acos((1+(x*x)-((d*d)*(1-x)*(1-x)))/(2*x));
  //return 0;
}

inline double PhononAbsorptionProcess::makeTDeviation(double d, double x){
  //change in T-phonon propagation direction after decay (L->L+T process)
  
  return acos((1-x*x+d*d*(1-x)*(1-x))/(2*d*(1-x)));
  //return 0;
}

inline double PhononAbsorptionProcess::makeTTDeviation(double d, double x){
  //change in T-phonon propagation direction after decay (L->T+T process)

  return acos((1-d*d*(1-x)*(1-x)+d*d*x*x)/(2*d*x));
}

void PhononAbsorptionProcess::makeTTSecondaries(const G4Track& aTrack){
  //Generate daughter phonons from L->T+T process
   
  //d is the velocity ratio vL/vT
  G4double d=1.6338;
  G4double upperBound=(1+(1/d))/2;
  G4double lowerBound=(1-(1/d))/2;

  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in first T phonon
  G4double x = d*(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
  G4double p = 1.5*G4UniformRand();
  while(!(p<getTTDecayProb(d, x))){
    x=d*(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
    p=1.5*G4UniformRand(); 
  }
  x=x/d;
  
  //using energy fraction x to calculate daughter phonon directions
  G4double theta1=makeTTDeviation(d, (x));
  G4double theta2=makeTTDeviation(d, (1-x));
  G4ThreeVector dir1=((PhononTrackInformation*)(aTrack.GetUserInformation()))->getK();
  G4ThreeVector dir2=dir1;
  G4ThreeVector ran = G4RandomDirection();
  
  //while(dir1.cross(ran)==0) ran=G4RandomDirection();
  //dir1 = dir1.rotate(dir1.cross(ran),theta1);
  //dir2 = dir2.rotate(dir2.cross(ran),-theta2);
  G4double ph=G4UniformRand()*2*PI;
  dir1 = dir1.rotate(dir1.orthogonal(),theta1).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-theta2).rotate(dir2,ph);

  aParticleChange.SetNumberOfSecondaries(2);
  G4double E=aTrack.GetKineticEnergy();
  G4Track* sec1;
  G4Track* sec2;
  
  G4double probST = Lattice->getSTDOS()/(Lattice->getSTDOS()+Lattice->getFTDOS());

 //First secondary:Make FT or ST phonon, probability density is funciton of eqn of state
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    sec1 = new G4Track(new G4DynamicParticle(TPhononSlow::PhononDefinition(),Lattice->mapKtoVDir(1,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }else{
    sec1 = new G4Track(new G4DynamicParticle(TPhononFast::PhononDefinition(),Lattice->mapKtoVDir(2,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }

 //Second secondary:Make FT or ST phonon, probability density is funciton of eqn of state
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    sec2 = new G4Track(new G4DynamicParticle(TPhononSlow::PhononDefinition(),Lattice->mapKtoVDir(1,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }else{
    sec2 = new G4Track(new G4DynamicParticle(TPhononFast::PhononDefinition(),Lattice->mapKtoVDir(2,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }

  //set the k-vectors for the two secondaries and add them to the process
  sec1->SetUserInformation(new PhononTrackInformation(dir1));
  sec2->SetUserInformation(new PhononTrackInformation(dir2));
  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2);


}

void PhononAbsorptionProcess::makeLTSecondaries(const G4Track& aTrack){

  //d is the velocity ratio vL/v
  G4double d=1.6338;
  G4double upperBound=1;
  G4double lowerBound=(d-1)/(d+1);
  
  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in L phonon
  G4double x =(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
  G4double p = 4.0*G4UniformRand();
  while(!(p<getLTDecayProb(d, x))){
    x=(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
    p=4.0*G4UniformRand(); //4.0 is about the max in the probability density function
  }

  //using energy fraction x to calculate daughter phonon directions
  G4double thetaL=makeLDeviation(d, x);
  G4double thetaT=makeTDeviation(d, x);;
  G4ThreeVector dir1=((PhononTrackInformation*)(aTrack.GetUserInformation()))->getK();
  G4ThreeVector dir2=dir1;

  G4double ph=G4UniformRand()*2*PI;
  dir1 = dir1.rotate(dir1.orthogonal(),thetaL).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-thetaT).rotate(dir2,ph);



  aParticleChange.SetNumberOfSecondaries(2);
  G4double E=aTrack.GetKineticEnergy();

  G4Track* sec1 = new G4Track(new G4DynamicParticle(LPhonon::PhononDefinition(),Lattice->mapKtoVDir(0,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );

  G4Track* sec2;

  //Make FT or ST phonon, probability density is funciton of eqn of state
  G4double probST = Lattice->getSTDOS()/(Lattice->getSTDOS()+Lattice->getFTDOS());
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    
    sec2 = new G4Track(new G4DynamicParticle(TPhononSlow::PhononDefinition(),Lattice->mapKtoVDir(1,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  
  }else{
    
    sec2 = new G4Track(new G4DynamicParticle(TPhononFast::PhononDefinition(),Lattice->mapKtoVDir(2,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  
  }

  sec1->SetUserInformation(new PhononTrackInformation(dir1));
  sec2->SetUserInformation(new PhononTrackInformation(dir2));
  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2);

    
}
