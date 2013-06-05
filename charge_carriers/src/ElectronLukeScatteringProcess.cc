#include "ElectronLukeScatteringProcess.hh"
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
#include <fstream>
#include <iostream>
#include "DriftingElectronTrackInformation.hh"


#include "G4Geantino.hh"


/*ElectronLukeScatteringProcess::ElectronLukeScatteringProcess(const G4String& aName)
  : G4VDiscreteProcess(aName)
{

  G4cout<<"\nWARNING:\n\tElectronLukeScatteringProcess::Constructor: Called constructor with pointer to step limitier! Process will not be properly initialised!!!"<<endl;
  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}
*/

ElectronLukeScatteringProcess::ElectronLukeScatteringProcess()
 :G4VDiscreteProcess("ElectronLukeScattering")
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

ElectronLukeScatteringProcess::~ElectronLukeScatteringProcess()
{ ; }

ElectronLukeScatteringProcess::ElectronLukeScatteringProcess(ElectronLukeScatteringProcess& right)
  : G4VDiscreteProcess(right)
{ ; }

G4double ElectronLukeScatteringProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  G4double me=9.1e-31*kg;
  G4double mx=0.081*me;
  G4double my=0.081*me;
  G4double mz=1.58*me;

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

  G4double velocity = aTrack.GetStep()->GetPostStepPoint()->GetVelocity();
  G4ThreeVector velDirXYZ = 
aTrack.GetStep()->GetPostStepPoint()->GetMomentumDirection().unit();
  G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
  ////////////////////////////////////
  // Debug: fix velocity (vector)
  //velocity = 107758.4471*m/s;
  //velDirXYZ = G4ThreeVector(-4.6293e-12,58986.1538,-90180.4668).unit();
  ///////////////////////////////////

  //  the following block computes velocity,kmag in HV space
  G4ThreeVector velDirLLT = normalToValley.TransformAxis(velDirXYZ).unit();
  G4ThreeVector velDirHV = G4ThreeVector(velDirLLT.getX()*(1/1.2172),
			 velDirLLT.getY()*(1/1.2172),
			 velDirLLT.getZ()*(1/0.27559)
			 );
  G4double konst = velDirHV.mag();
  G4double kmag = (velocity*massHole / hbar) * konst; 
  //G4cout <<  "Dan's kmag = " <<  kmag*m <<  G4endl;
  //G4cout <<  "Dan's konst = " <<  konst <<  G4endl;
  
  
  //G4ThreeVector k = 
//postStepPoint->GetMomentum()/hbar_Planck;
  /*G4ThreeVector k = 
postStepPoint->GetVelocity()*postStepPoint->GetMass()/c_squared/hbar_Planck*
  postStepPoint->GetMomentumDirection().unit();
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV = G4ThreeVector(k_valley[0]/1.2172, k_valley[1]/1.2172, 
				      k_valley[2]/0.27559);
				      
  G4double kmag = k_HV.mag();
  */
  

  *condition = NotForced;
  //G4cout <<  "k = " <<  k.mag()*m <<  G4endl;
  //G4cout << "kmag = " << kmag*m << G4endl;
  //G4cout << "ksound_Hole = " << ksound_Hole*m << G4endl;

  if (kmag<ksound_Hole) return DBL_MAX;  

  G4double tau =  1 / (
		       velLong / (3*l0Hole)
		       * (kmag / ksound_Hole)*(kmag/ksound_Hole)
		       * ((1- ksound_Hole /kmag))
		       * ((1- ksound_Hole /kmag))
		       * ((1- ksound_Hole /kmag))
		       );
  
  //////////////////////////
  // Debug
  //tau = 0.01*ns;
  /////////////////////////
  
  G4double mfp = tau*velocity;
  /*
  G4cout<<"\nElectronLukeScatteringProcess::GetMeanFreePath:";
  G4cout<<"\n\tvelocity vector: \t"<<velDirXYZ.unit();
  G4cout<<"\n\tvelocity XYZ:\t"<<velocity / (m/s)<<" m/s";
  G4cout<<"\n\tvelocity HV:\t"<<velocity*konst / (m/s)<<" m/s";
  G4cout<<"\n\tkmag:\t\t:"<<kmag;
  G4cout<<"\n\tvelLong:\t"<<velLong /(m/s)<<" m/s";
  G4cout<<"\n\tl0Hole:\t"<<l0Hole/m<<" m";
  G4cout<<"\n\tvelLong/(3l0Hole):\t"<<velLong / (3*l0Hole) /(1/s)<<" 1/s";
  G4cout<<"\n\tkmag / ksound_Hole:\t"<<kmag / ksound_Hole;
  G4cout<<"\n\t(1- ksound_Hole /kmag):\t"<<(1- ksound_Hole /kmag);
  G4cout<<"\n\ttau:\t\t"<<tau / (s) <<" (s)";
  */
  //G4cout<<"\n\tmfp:\t\t"<<mfp / mm<<" mm";
  //G4cout<<endl;
  
  return mfp;
 
}

G4VParticleChange* ElectronLukeScatteringProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  /*
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  G4RotationMatrix trix;
  int valley = ((DriftingElectronTrackInformation*) 
aTrack.GetUserInformation())->getValley();
  switch(valley){
      case 1:
          trix = G4RotationMatrix(PI/4, PI/4, PI/4);
          break;
      case 2:
          trix = G4RotationMatrix(-PI/4, PI/4, PI/4);
          break;
      case 3:
          trix = G4RotationMatrix(PI/4, -PI/4, PI/4);
          break;
      case 4:
          trix = G4RotationMatrix(-PI/4, -PI/4, PI/4);
          break;
  }
  
  normalToValley= G4AffineTransform(trix);
  valleyToNormal= G4AffineTransform(trix).Inverse();
  
  G4ThreeVector k = postStepPoint->GetMomentum()/postStepPoint->GetGamma()/hbar;
  //G4cout <<  "Kinetic Energy: " <<  postStepPoint->GetKineticEnergy()/eV <<  
//G4endl;
//G4cout <<  "Total Energy: " <<  postStepPoint->GetTotalEnergy()/eV << G4endl;
  //G4cout <<  "k^2 h^2/2m: " <<  
//k.mag2()*hbar*hbar/2/postStepPoint->GetMass()/eV <<  G4endl;
  G4ThreeVector k_valley = normalToValley.TransformPoint(k);
  G4ThreeVector k_HV = G4ThreeVector(k_valley[0]/1.2172, k_valley[1]/1.2172, 
				      k_valley[2]/0.27559);
				      
  G4double kmag = k_HV.mag();
  G4double theta_phonon = MakeTheta(kmag, ksound_Hole);

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


  G4double qEnergy = (
		      4 * hbar*hbar*ksound_Hole*ksound_Hole / 2 
		      /postStepPoint->GetMass() * (kmag / ksound_Hole 
				     * cos(theta_phonon) - 1 
				     )
		      );
      
  G4ThreeVector k_HV_new = k_HV.rotate(k_HV.orthogonal(), theta_charge);
  G4double phi_charge =  G4UniformRand()*2*pi;
  k_HV_new = k_HV_new.rotate(k_HV_new.unit(), phi_charge);
  G4ThreeVector k_valley_new = G4ThreeVector(k_HV_new[0]*1.2172, 
					     k_HV_new[1]*1.2172, 
					     k_HV_new[2]*0.27559);
  G4ThreeVector k_new = valleyToNormal.TransformPoint(k_valley_new);
  G4double E_new = (k_HV_new.mag2()*hbar*hbar)
		  /2/postStepPoint->GetMass()- qEnergy;
G4cout << qEnergy <<  G4endl;
G4cout <<  E_new <<  G4endl;
  aParticleChange.ProposeMomentumDirection(k_new.unit());
  aParticleChange.ProposeEnergy(E_new);
  ResetNumberOfInteractionLengthLeft();    

  return &aParticleChange;
  */

  aParticleChange.Initialize(aTrack);  
  //Do nothing other than re-calculate mfp when step limit reached or leaving 
  //volume
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  //if(  
   //  (postStepPoint->GetProcessDefinedStep()==stepLimiter)
  //   ||(postStepPoint->GetStepStatus()==fGeomBoundary)
  //   )
  //  {     
  //     return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
 //   }  
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
  
  //The following block calculates the scaling constant between HV and real space
  
  G4ThreeVector DirXYZ = aTrack.GetMomentumDirection();  
  G4double velocityXYZ = aStep.GetPostStepPoint()->GetVelocity();
  ////////////////////////////////
  // Debug: fix momentum vector and velocity coming in
  //DirXYZ = G4ThreeVector(-4.771e-17  ,   0.54739  ,  -0.83688);//(-4.771e-17 ,    0.54739  ,  -0.83688);
  //velocityXYZ = 44083.0011 * m/s;
  ///////////////////////////////

  G4ThreeVector DirLLT = normalToValley.TransformAxis(DirXYZ).unit();
  G4ThreeVector DirHV = G4ThreeVector(DirLLT.getX()*(1/1.2172),
		      DirLLT.getY()*(1/1.2172),
		      DirLLT.getZ()*(1/0.27559)
		      );
  G4double konst = DirHV.mag();
  //  DirHV = valleyToNormal.TransformAxis(DirHV).unit();
  
  G4double velocityHV = velocityXYZ * konst;
  G4double kmag = velocityXYZ*massHole / hbar * konst;
  //G4cout <<  "konst: " <<  konst <<  G4endl;
  //G4cout <<  "velocityXYZ: " <<  velocityXYZ<<  G4endl;

  G4double theta_phonon= MakeTheta(kmag, ksound_Hole);

  /////////////////////////////////////////
  //DEBUG - FIXING theta_phonon FOR TESTING
  //theta_phonon = 1.1543;//pi/2;//1.1543;
  //kmag = 3.7949e7*(1/m);
  ////////////////////////////////////////

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
      //G4cout <<  "ksound_Hole " <<  ksound_Hole <<  G4endl;
      //G4cout <<  "massHole " <<  massHole <<  G4endl;
      //G4cout <<  "kmag " <<  kmag <<  G4endl;
      //G4cout <<  "theta_phonon " <<  theta_phonon <<  G4endl;
 
  G4ThreeVector orth = DirHV.orthogonal();
  G4ThreeVector DirHV_new = G4ThreeVector(DirHV).rotate(orth, theta_charge).unit();
  G4double phi_charge =  G4UniformRand()*2*pi;
  ///////////////////////////////
  // Debug: Set phi = 0
  //phi_charge = 4.77823;
  //////////////////////////////
  DirHV_new.rotate(DirHV, phi_charge);
  G4ThreeVector DirXYZ_new = G4ThreeVector(DirHV_new.getX()*1.2172,
					   DirHV_new.getY()*1.2172,
					   DirHV_new.getZ()*0.27559
					   );
  konst = DirXYZ_new.mag();
 // newDir = G4ThreeVector(newDir.getX()*1/0.2756,
//		         newDir.getY()*1/1.2172,
//			 newDir.getZ()*1/1.2172
//			 );
  DirXYZ_new = valleyToNormal.TransformAxis(DirXYZ_new).unit();
  
  //konst = 1.000*valleyToNormal.TransformAxis(newDir).mag();
  //G4cout<<"\nElectronLukeScatteringProcess::PostStepDoIt: new const: "<<konst;

  G4double kmag_new = sqrt(kmag*kmag
				  -(qEnergy*2*massHole/hbar/hbar)
			   );
   G4double v_new = hbar/massHole * konst*kmag_new;
   G4double E_new = 0.5*massHole*v_new*v_new;
   G4double E_old = 0.5*massHole*velocityXYZ*velocityXYZ;
   //G4cout <<  "kmag_new " <<  kmag_new <<  G4endl;
   //G4cout <<  "kmag " <<  kmag<<  G4endl;
   //G4cout <<  "qEnergy " <<  qEnergy <<  G4endl;
   //G4cout <<  "v_new " <<  v_new <<  G4endl;
   //G4cout <<  "E_new " <<  E_new <<  G4endl;
   //G4cout <<  "E_old " <<  E_old <<  G4endl;
   //G4cout<<"\n\nElectronLukeScattering::PostStepDoIt:\n\tPrior to 
//scattering:";
   //G4cout<<"\n\tvelocity vector XYZ:\t"<<DirXYZ.unit();
   //G4cout<<"\n\tvelocity vector LLT:\t"<<DirLLT.unit();
   //G4cout<<"\n\tvelocity vector HV:\t"<<DirHV.unit();
   //G4cout<<"\n\tvelocity in XYZ:\t"<<velocityXYZ/(m/s)<<"m/s";
   //G4cout<<"\n\tvelocity in HV:\t"<<velocityHV/(m/s)<<"m/s";
   //G4cout<<"\n\tkmag: \t\t"<<kmag / (1/m)<<" [1/m]";
   //G4cout<<"\nPostScattering:";
   //G4cout<<"\n\tnew velocity vector XYZ:\t"<<DirXYZ_new.unit();
   //G4cout<<"\n\tnew velocity vector HV:\t"<<DirHV_new.unit();
   //G4cout<<"\n\tnew kmag: \t\t"<<kmag_new / (1/m)<<" [1/m]";
   //G4cout<<"\n\tnew velocity in HV:\t"<<kmag_new*hbar/massHole/(m/s)<<"m/s";
   //G4cout<<"\n\tnew velocity in XYZ:\t"<<v_new/(m/s)<<"m/s";
   //G4cout<<"\n\talternative new velocity in XYZ:\t"<< 
//kmag_new*konst*hbar/massHole/(m/s) /(m/s)<<"m/s";
   //G4cout<<"\n\tqEnergy:\t"<<qEnergy/eV<<" [eV]";
   //G4cout<<"\n\ttheta phonon:\t\t"<<theta_phonon/rad<<" [rad]";
   //G4cout<<"\n\ttheta charge:\t\t"<<theta_charge/rad<<" [rad]";
   //G4cout<<"\n\tphi charge:\t\t"<<phi_charge/rad<<" [rad]";
   //G4cout<<"\n\told Ekin:\t"<<aTrack.GetKineticEnergy()/eV<<" [eV]";
   //G4cout<<"\n\told Ekin from velocity:\t"<<E_old/eV<<" [eV]";
   //G4cout<<"\n\tnew Ekin:\t"<<E_new /eV<<" [eV]";
   //G4cout<<endl;
   
   
   //G4cout<<"\n\tnew velocity vector HV:\t"<<DirHV_new.unit()<<" phi: "<<phi_charge/rad;
  //G4cout<<"\n\tkonst:\t"<<konst;
  //G4cout<<"\n\tkmag: \t\t"<<kmag / (1/m)<<" [1/m]";
  //G4cout<<"\n\tangle:\t\t"<<theta_phonon/rad<<" [rad]";
  //G4cout<<"\n\tqenergy:\t"<<qEnergy/eV<<" [eV]";
  //G4cout<<"\n\tkmag_new:\t"<<kmag_new / (1/m)<<" [1/m]";
  //G4cout<<"\n\tKE_old:\t"<<hbar*hbar*kmag*kmag / (2*massHole) /eV<<" [eV]";
  //G4cout<<"\n\tKE_new:\t"<<hbar*hbar*kmag_new*kmag_new / (2*massHole) /eV<<" 
//[eV]";
  //G4cout<<"\n\tv_new:\t"<<v_new/(m/s)<<" [m/s]";
  //G4cout<<"\n\tKE_new2:\t"<<E_new /eV<<" [eV]";
   


  G4double kineticEnergy =( E_new);//hbar*hbar*kmag_new*kmag_new / (2*massHole) );
  //G4cout<<"\nElectronLukeScatteringProcess::PostStepDoIt: kineticEnergy: "<<kineticEnergy/eV;
  //if(kineticEnergy<=0.0) kineticEnergy=DBL_MIN;

  aParticleChange.ProposeMomentumDirection(DirXYZ_new.unit());
  //aParticleChange.ProposeMomentumDirection(G4ThreeVector(1,0,0));
  /////////////////
  // Debug: do not change direction of propagation
  //aParticleChange.ProposeMomentumDirection(G4RandomDirection());
  ////////////////
  //aParticleChange.ProposeEnergy(1e-12*eV);
  aParticleChange.ProposeEnergy(E_new);
  ResetNumberOfInteractionLengthLeft();    

  //G4ThreeVector newPosn = aTrack.GetPosition();
  //std::ofstream epositions;
  //epositions.open("e-positions.txt", std::ofstream::app);

  //epositions << newPosn.getX() << " " << newPosn.getY() << " " << newPosn.getZ() << "\n";
  //epositions.close();
  return &aParticleChange;
}

G4double ElectronLukeScatteringProcess::MakeTheta(G4double& k, G4double& ks){
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

  bool first=false;
    
  while(pValue>pDensity){
    theta=G4UniformRand()*thetaMax;
    pValue=G4UniformRand();// *pValMax;//  *(1+2*ks/k+(ks/k)*(ks/k));   
    pDensity = (cos(theta)-(ks/k))*(cos(theta)-(ks/k))*sin(theta);//need to multiply by unit 's' to make it dimensionless
    if(pDensity>pValMax) G4cout<<"\nLukeScattering::PostStepDoIt: Error: pDensity should never exceed pValMax "<<pValMax;
    
    
    //    if(!first){
    //G4cout<<"\nLukeScatteringProcess::MakeTheta: pDensity calculated as: "<<pDensity;
    // G4cout<<"\n\tpValue: "<<pValue;
    // G4cout<<"\n\ttheta: "<<theta/rad;
    // G4cout<<"\n\tk: "<<k*m;
    // G4cout<<"\n\tks: "<<ks*m;
    // G4cout<<endl;
    //}
    //first=1;
    
    }

  return theta;
}

G4double ElectronLukeScatteringProcess::MakePhi(G4double& k,G4double& ks, G4double& theta){


  G4double phi=acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));

  return phi;
}

G4bool ElectronLukeScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  
return((&aPD==DriftingElectron::Definition()));
}
