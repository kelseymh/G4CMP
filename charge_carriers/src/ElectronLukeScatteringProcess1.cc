#include "ElectronLukeScatteringProcess1.hh"
#include "DriftingElectron.hh"
#include "DriftingHole.hh"

#include "LPhonon.hh"
#include "PhononTrackInformation.hh"
#include "LatticeManager2.hh"
#include "PhysicalLattice.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4TransportationManager.hh"
#include <fstream>
#include <iostream>
#include "DriftingElectronTrackInformation.hh"
#include "G4Field.hh"

#include "G4FieldManager.hh"

#include "G4Geantino.hh"


/*ElectronLukeScatteringProcess::ElectronLukeScatteringProcess(const G4String& aName)
  : G4VProcess(aName)
{

  G4cout<<"\nWARNING:\n\tElectronLukeScatteringProcess::Constructor: Called constructor with pointer to step limitier! Process will not be properly initialised!!!"<<endl;
  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}
*/

ElectronLukeScatteringProcess1::ElectronLukeScatteringProcess1()
 :G4VProcess("ElectronLukeScattering")
{
  velLong=5324.2077*m/s;
  l0Hole = 257e-6*m;
  massFreeElectron=9.1e-31*kg;
  massHole=3/(1/1.58+2/.081)*massFreeElectron;
  hbar= 6.5821e-16*eV*s;
  //ksound_Hole=5512124.0235/m; //velLong*massHole / hbar;//1.6306e7/m;
  ksound_Hole = velLong*massHole/hbar_Planck;
  G4cout<<"\nElectronLukeScatteringProcess::Constructor:ksound_Hole = "
	<<ksound_Hole*m<<" /m"<<G4endl;

  /*
    normalToValley= G4AffineTransform
    (G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532));
  valleyToNormal= G4AffineTransform
    (G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532)).Inverse();
  */

  /*
  HVMatrix = G4RotationMatrix(G4ThreeVector(1.2172, 0     ,     0),
				      G4ThreeVector(0     , 1.2172,     0),
				      G4ThreeVector(0     , 0     ,0.2756)
				      );
  HVTransform = G4AffineTransform(HVMatrix);
  */

  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}

ElectronLukeScatteringProcess1::~ElectronLukeScatteringProcess1()
{ ; }

ElectronLukeScatteringProcess1::ElectronLukeScatteringProcess1(
ElectronLukeScatteringProcess1& right)
  : G4VProcess(right)
{ ; }

G4double 
ElectronLukeScatteringProcess1::PostStepGetPhysicalInteractionLength(const 
G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  G4double me=9.1e-31*kg;
  G4double mx=0.081*me;
  G4double my=0.081*me;
  G4double mz=1.58*me;
  G4double mc=.12*me;

  G4RotationMatrix trix;
  int valley = ((DriftingElectronTrackInformation*) aTrack.GetUserInformation())->getValley();
  switch(valley){
      case 1:
          trix = G4RotationMatrix(-PI/4, -PI/4, PI/4);
          break;
      case 2:
          trix = G4RotationMatrix(PI/4, -PI/4, -PI/4);
          break;
      case 3:
          trix = G4RotationMatrix(-PI/4, PI/4, PI/4);
          break;
      case 4:
          trix = G4RotationMatrix(PI/4, PI/4, -PI/4);
          break;
  }
  //trix = G4RotationMatrix(G4ThreeVector(0   ,   0.8165   ,   0.57735),
   //					   G4ThreeVector( 1 ,  0    ,        0),
	//				   G4ThreeVector(0   ,   0.57735 ,     -0.8165)
//					   );
  
  normalToValley= G4AffineTransform(trix);
  valleyToNormal= G4AffineTransform(trix).Inverse();

  
  G4ThreeVector k = aTrack.GetMomentum()/hbar_Planck;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV = G4ThreeVector(k_valley[0]*1.2172, k_valley[1]*1.2172, 
				      k_valley[2]*0.27559);
  //G4double kmag = k_HV.mag();
    G4FieldManager* fMan 
= G4TransportationManager::GetTransportationManager()->GetFieldManager();
    const G4Field* field = fMan->GetDetectorField();


    G4ThreeVector posVec = aTrack.GetPosition();
    G4double position[4] = {posVec[0],posVec[1],posVec[2],0};
    G4double fieldVal[6];

    field->GetFieldValue(position,fieldVal);
    G4ThreeVector Efield = G4ThreeVector(fieldVal[4], fieldVal[5], fieldVal[6]);
    G4ThreeVector Efield_valley = normalToValley.TransformPoint(Efield);
    G4ThreeVector Efield_HV = G4ThreeVector( Efield_valley[0]*1.2172, 
					     Efield_valley[1]*1.2172, 
					     Efield_valley[2]*0.27559);
    G4double Emag = Efield_HV.mag();
    G4ThreeVector v = G4ThreeVector(k_HV[0]*1.2172, k_HV[1]*1.2172, 
k_HV[2]*.27559)*hbar_Planck/mc;
    G4ThreeVector v_HV = k_HV*hbar_Planck/mc;
    G4double vmag[3] = {v_HV[0], v_HV[1], v_HV[2]};
    G4double e = electron_charge;
    G4double kp[8];
    G4double kmag[3] = {k_HV[0], k_HV[1], k_HV[2]};
    G4double ks[8];
    ks[0] = ksound_Hole;
    for(int i=1;i<8;++i) ks[i]=ks[i-1]*ksound_Hole;
    G4double v2;
    G4double vl2 = velLong*velLong;
    G4double v4;
    G4double v6;
    G4double E[3] = {Efield_HV[0], Efield_HV[1], Efield_HV[2]};
    G4ThreeVector step;
    
    /*
    for (int i=0;i<3;++i)
    {
    kp[0] = kmag[i];
    for(int j=1;j<8;++j) kp[j]=kp[j-1]*kmag[i];
    v2 = vmag[i]*vmag[i];
    v4 = v2*v2;
    v6 = v4*v2;
    
    G4double denom = 
e*e*E[i]*E[i]*vmag[i]*(hbar_Planck*kp[3]*(kp[2]-ks[0]*(3*kp[1]+ks[0]*(-3*kp[0]+
ks [ 0
] ) ) ) - 2*hbar_Planck*(kp[4]-3*kp[2]*ks[1]+2*kp[1]*ks[2])*mc*v2 + 
(kp[2]-9*kp[0]*ks[1]+8*ks[2])*mc*mc*v4);

    G4double t1 = 3*hbar_Planck*hbar_Planck*kp[3]*
(kp[5]-20*kp[2]*ks[2]+15*kp[1]*ks[3]-6*kp[0]*ks[4]+ks[5]+3*kp[3]*ks 
[0]*(-2*kp[0]+5*ks[0]));
    G4double t2 = 
-2*hbar_Planck*(kp[7]-3*kp[6]*ks[0]+5*kp[3]*(2*kp[0]-3*ks[0])*ks[2]+9*kp[2]*ks[4
]-2*k[1]*ks[5])*mc*v2;
    G4double t3 = 
3*kp[5]+24*kp[2]*ks[2]-31*kp[1]*ks[3]+18*kp[0]*ks[4]-4*ks[5]-2*kp[3]*ks[0]*(kp[0
]+3*ks[0] );
    G4double t4 = hbar_Planck*(kp[3]*(kp[0]-3*ks[0])+3*kp[2]*ks[1]-kp[1]*ks[2]) 
+ (kp[0]-ks[0])*(kp[0]-ks[0])*(kp[0]+2*ks[0])*mc*v2;
    G4double t5 = 
4*hbar_Planck*hbar_Planck*kp[7]*ks[1]*l0Hole*l0Hole*mc*mc*v2*velLong*velLong;
    G4double t0 = 
hbar_Planck*mc*(electron_charge*E[i]*v2*vmag[i]*(hbar_Planck*kp[3]*(kp[2]-ks[0]*
( 3 * kp[1]-3*kp[0]*ks[0]+ks[1])) + (kp[4]-3*kp[2]*ks[1]+2*kp[1]*ks[2])*mc*v2) 
- 2*hbar_Planck*kp[5]*ks[0]*l0Hole*mc*v4*velLong);
    
    step[i] = 
(t0+sqrt(hbar_Planck*hbar_Planck*kp[3]*mc*mc*v2*(electron_charge*E[i]*E[i]*(t1+
t2 +t3)*mc*mc*v4- 
4*electron_charge*E[i]*hbar_Planck*kp[3]*ks[0]*l0Hole*mc*vmag[i]*t4*velLong+t5 
) ))/denom;

if (!(step[i] > 0))
     step[i] = 
(t0+sqrt(hbar_Planck*hbar_Planck*kp[3]*mc*mc*v2*(electron_charge*E[i]*E[i]*(t1+
t2 +t3)*mc*mc*v4- 
4*electron_charge*E[i]*hbar_Planck*kp[3]*ks[0]*l0Hole*mc*vmag[i]*t4*velLong+t5 
) ))/denom;
    }
    */
    
  
  *condition = NotForced;
  //G4cout << "Luke recalc mfp" <<  G4endl;
  
  G4cout << step.mag()/m << G4endl;
  return step.mag();
 
}

G4VParticleChange* ElectronLukeScatteringProcess1::PostStepDoIt(const G4Track& 
aTrack, const G4Step& aStep)
{
 G4cout << "Scatter" << G4endl;
  G4double me=9.1e-31*kg;
  G4double mx=0.081*me;
  G4double my=0.081*me;
  G4double mz=1.58*me;
  G4double mc = .12*me;
  aParticleChange.Initialize(aTrack); 
  G4RotationMatrix trix;
  int valley = ((DriftingElectronTrackInformation*) 
aTrack.GetUserInformation())->getValley();
  switch(valley){
      case 1:
          trix = G4RotationMatrix(-PI/4, -PI/4, PI/4);
          break;
      case 2:
          trix = G4RotationMatrix(PI/4, -PI/4, -PI/4);
          break;
      case 3:
          trix = G4RotationMatrix(-PI/4, PI/4, PI/4);
          break;
      case 4:
          trix = G4RotationMatrix(PI/4, PI/4, -PI/4);
          break;
  }
  
  normalToValley= G4AffineTransform(trix);
  valleyToNormal= G4AffineTransform(trix).Inverse();
  
  G4ThreeVector v = aTrack.GetVelocity()*aTrack.GetMomentumDirection();
  G4ThreeVector v_Valley = normalToValley.TransformPoint(v);
  G4ThreeVector v_HV = G4ThreeVector(v_Valley[0]/1.217, 
					v_Valley[1]/1.217, 
					v_Valley[2]/0.27559);
  G4ThreeVector k = aTrack.GetMomentum()/hbar_Planck/c_light;
  G4ThreeVector k_Valley = normalToValley.TransformPoint(k);
/*  G4cout <<  (hbar_Planck*hbar_Planck*k_Valley[0]*k_Valley[0]/2/mx +
	    hbar_Planck*hbar_Planck*k_Valley[1]*k_Valley[1]/2/my +
	    hbar_Planck*hbar_Planck*k_Valley[2]*k_Valley[2]/2/mz)/MeV <<  
G4endl;
*/
//G4cout <<  "Kinetic Energy = " << aTrack.GetKineticEnergy()/MeV <<  G4endl;
	    
  G4ThreeVector k_HV = G4ThreeVector(k_Valley[0]*1.217, 
					k_Valley[1]*1.217, 
					k_Valley[2]*0.27559);
//G4ThreeVector k_HV = v_HV*mc/hbar_Planck;
   // k_HV.setMag(k.mag());
  G4double kmag = k_HV.mag(); 
  //Do nothing other than re-calculate mfp when step limit reached or leaving 
  //volume
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if(  
     (postStepPoint->GetProcessDefinedStep()==stepLimiter)
     ||(postStepPoint->GetStepStatus()==fGeomBoundary
       || kmag <= ksound_Hole)
     )
    {     
       return &aParticleChange;
    }  
    G4cout <<  "Scatter" <<  G4endl;

  G4double theta_phonon= MakeTheta(kmag, ksound_Hole);

  G4double theta_charge=acos( 
			     (kmag*kmag - 2*ksound_Hole
			      *(kmag*cos(theta_phonon) - ksound_Hole) 
			      - 2 * (kmag*cos(theta_phonon) - ksound_Hole)
			      * (kmag*cos(theta_phonon) - ksound_Hole)
			      ) 
			     / kmag
			     / (sqrt(
				     kmag*kmag - 4*ksound_Hole
				     *(kmag*cos(theta_phonon) - ksound_Hole)
				     )) 
			     );


  G4double qEnergy = (4 * hbar*hbar *ksound_Hole*ksound_Hole/2/massHole  
		      * (kmag / ksound_Hole * cos(theta_phonon) - 1) );
	//G4cout <<  "qEnergy = " << qEnergy/MeV <<  G4endl;
 
 G4ThreeVector kdir = k_HV.unit();
 G4ThreeVector kdir_new = kdir.rotate(kdir.orthogonal(), theta_charge);
  G4double phi_charge =  G4UniformRand()*2*pi;
  kdir_new.rotate(kdir, phi_charge);
  G4double q = sqrt(qEnergy*2*massHole/hbar_Planck/hbar_Planck);
  k_HV = 
kdir_new*(sqrt(k_HV.mag2()+q*q-q*k_HV.mag()*cos(theta_phonon)));
  
  k_Valley = G4ThreeVector(k_HV[0]/1.217, 
		    k_HV[1]/1.217, 
		    k_HV[2]/0.27559);
  k = valleyToNormal.TransformPoint(k_Valley);
  G4ThreeVector v_new = hbar_Planck/mc*G4ThreeVector(k_HV[0]*1.217, 
k_HV[1]*1.217, k_HV[2]*0.27559);
  G4double E_new = hbar_Planck*hbar_Planck*k_Valley[0]*k_Valley[0]/2/mx 
		 + hbar_Planck*hbar_Planck*k_Valley[1]*k_Valley[1]/2/my
		 + hbar_Planck*hbar_Planck*k_Valley[2]*k_Valley[2]/2/mz;
  //G4double E_new = .5*v.mag2()*me;
    //G4cout << aTrack.GetKineticEnergy()/MeV <<  G4endl;
    //G4cout << "Enew = " <<  E_new/MeV <<  G4endl;

    G4cout << v_new.mag()*s/m<< G4endl;
  aParticleChange.ProposeMomentumDirection(v_new.unit());
  //aParticleChange.ProposeEnergy(E_new);
  aParticleChange.ProposeVelocity(v_new.mag());
  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}

G4double ElectronLukeScatteringProcess1::MakeTheta(G4double& k, G4double& ks){
  //G4double u = G4UniformRand();
 
  //Analytical method to compute theta
  
 /* 

  double base = -(u-1)+3*(u-1)*(ks/k)-3*(u-1)*(ks/k)*(ks/k)+(u-1)*(ks/k)*(ks/k)*(ks/k);
  double exponent =1.0/3.0;
  G4double operand = ks/k+pow(base, exponent);   

  if(operand>1.0) {
    // G4cout<<"\nTruncating operand from"<<operand<<" to 1.0";
    operand=1.0;
  }
  if(operand<0.0) G4cout<<"\noperand out of range: operand = "<<operand;
  G4double theta=acos(operand);
  //  G4cout<<"\n"<<theta;

  if(acos(ks/k)<theta) G4cout<<"\n  THETA OUT OF BOUNDS!!! (theta>acos(ks/k)) theta:"<<theta<<" acos(ks/k):"<<acos(ks/k);
  if(PI/2<=theta) G4cout<<"\n THETA OUT OF BOUNDS!!! (pi/2 < theta)";

  */
  //Rejection method for determining theta
  G4double theta;
  G4double pValue=2;
  G4double pDensity = 1;

  G4double thetaMax;
  G4double pValMax;
  thetaMax=acos(ks/k);
  pValMax=(1-(ks/k))*(1-(ks/k));
  bool first = true;

  while(pValue>pDensity){
    pValue=G4UniformRand();//*pValMax;//  *(1+2*ks/k+(ks/k)*(ks/k));   
    theta=G4UniformRand()*thetaMax;
    pDensity = (cos(theta)-(ks/k))*(cos(theta)-(ks/k))*sin(theta);//need to multiply by unit 's' to make it dimensionless
    if(pDensity>pValMax) G4cout<<"\nLukeScattering::PostStepDoIt: Error: pDensity should never exceed pValMax "<<pValMax;
    
    
    //    if(!first){
    //G4cout<<"\nLukeScatteringProcess::MakeTheta: pDensity calculated as: 
//"<<pDensity;
  //   G4cout<<"\n\tpValue: "<<pValue;
    // G4cout<<"\n\ttheta: "<<theta/rad;
    // G4cout<<"\n\tk: "<<k*m;
   //  G4cout<<"\n\tks: "<<ks*m;
   //  G4cout<<endl;
    //}
    first=false;
    
  }
 /*
    //Analytical Solution for theta
    // Let G(y) = F^(-1)(x)
    // If u(1)...u(n) are randoms from the uniform dist. (0, 1)
    // then G(u(1))...G(u(n)) is from the cdf F(x)
    G4double theta;
    G4double u;
    G4double vs = velLong;
    //G4cout <<  "vs = " <<  vs <<  G4endl;
    G4double l0 = l0Hole;
    //G4cout <<  "l0 = " <<  l0 <<  G4endl;
    G4double theta_max = acos(ks/k);
    //G4cout <<  "theta_max = " <<  theta_max <<  G4endl;
    //G4cout <<  "k = " <<  k*m <<  G4endl;
    //G4cout <<  "ks = " <<  ks*m <<  G4endl;
    if (theta_max > .01) //short-circuit if max angle is ~< 1 degree
    do {
    u = G4UniformRand();
    //G4cout <<  "u = " <<  u <<  G4endl;
    theta = acos(ks/k + 1/k/k/vs*pow(-3*u*(k*k*k*k)*(ks*ks)*l0*(vs*vs) + 
		(k*k*k*k*k*k)*(vs*vs*vs) - 3*(k*k*k*k*k)*ks*(vs*vs*vs) + 
		3*(k*k*k*k)*(ks*ks)*(vs*vs*vs) - (k*k*k)*(ks*ks*ks)*(vs*vs*vs), 
		1/3) );
    //theta = acos((pow(-(u-1)*k*k*k*k*k*k*(k-ks)*(k-ks)*(k-ks), 1/3) + 
                //k*k*ks)/k/k/k);
    //theta = acos(1/k/k*(k*ks+pow(-k*k*k*((u-1)*k*k*k - 3*(u-1)*k*k*ks + 
                //3*(u-1)*k*ks*ks + ks*ks*ks), 1/3)));
	//	G4cout <<  theta <<  G4endl;
    } while (theta>theta_max || isnan(theta) );
    else
        theta = 0;
	*/
  return theta;
}

G4double ElectronLukeScatteringProcess1::MakePhi(G4double& k,G4double& ks, 
G4double& theta){


  G4double phi=acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));

  return phi;
}

G4bool ElectronLukeScatteringProcess1::IsApplicable(const G4ParticleDefinition& 
aPD)
{
  
return((&aPD==DriftingElectron::Definition()));
}
