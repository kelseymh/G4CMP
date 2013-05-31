#include "LukeScatteringProcess.hh"
#include "DriftingElectron.hh"
#include "DriftingHole.hh"

#include "LPhonon.hh"
#include "PhononTrackInformation.hh"
#include "LatticeManager2.hh"
#include "PhysicalLattice.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"

#include "G4Geantino.hh"

#include "math.h"


/*LukeScatteringProcess::LukeScatteringProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{
if(verboseLevel>1){
G4cout<<GetProcessName()<<" is created "<<G4endl;
}
}
 */
LukeScatteringProcess::LukeScatteringProcess()
:G4VDiscreteProcess("LukeScattering")
{
    velLong=5400*m/s;
    l0Hole = 108e-6*m;
    massFreeElectron=9.1e-31*kg;
    massHole=0.35*massFreeElectron;
    ksound_Hole=1.6306e7/m;
    hbar= 6.5821e-16*eV*s;

    if(verboseLevel>1){
	G4cout<<GetProcessName()<<" is created "<<G4endl;
      }
}

LukeScatteringProcess::~LukeScatteringProcess()
{ ; }

LukeScatteringProcess::LukeScatteringProcess(LukeScatteringProcess& right)
: G4VDiscreteProcess(right)
{ ; }

G4double LukeScatteringProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
    G4StepPoint* stepPoint  = aTrack.GetStep()->GetPostStepPoint();
    G4double velocity = stepPoint->GetVelocity();
    G4double kmag = velocity*massHole / hbar;

    *condition = NotForced;

    if (kmag<=ksound_Hole)
    {
	//G4cout <<  "LukeScatteringProcess: DBL_MAX" <<  G4endl;
	return DBL_MAX;
    }

    G4double mfp= velocity / ( velLong / (3*l0Hole)
			    * (kmag / ksound_Hole)*(kmag / ksound_Hole)
			    * ((1- ksound_Hole /kmag))
			    * ((1- ksound_Hole /kmag))
			    * ((1- ksound_Hole /kmag)));


    //G4cout <<  "LukeScatteringProcess: " <<  mfp <<  G4endl;
    return mfp;

}

G4VParticleChange* LukeScatteringProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{

    aParticleChange.Initialize(aTrack);  

    //Do nothing other than re-calculate mfp when step limit reached or leaving volume
    G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
    //G4cout<<"\nLukeScattering::PostStepDoIt: Step limited by process " <<postStepPoint->GetProcessDefinedStep()->GetProcessName() <<  G4endl;


    /*if(	(postStepPoint->GetProcessDefinedStep()==stepLimiter)
	|| (postStepPoint->GetStepStatus()==fGeomBoundary))
	return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    */


    G4double velocity = aStep.GetPostStepPoint()->GetVelocity();
    G4double kmag = velocity*massHole / hbar;
    G4double theta_phonon=MakeTheta(kmag, ksound_Hole);
    G4double theta_charge=acos( 
					    (kmag*kmag - 2*ksound_Hole
					    *(kmag*cos(theta_phonon) - ksound_Hole) 
					    - 2 * (kmag*cos(theta_phonon) - ksound_Hole)
					    * (kmag*cos(theta_phonon) - ksound_Hole)) 
					    / kmag/ (sqrt(kmag*kmag - 4*ksound_Hole
					    *(kmag*cos(theta_phonon) - ksound_Hole))));

    G4double qEnergy = (4 * hbar*hbar *ksound_Hole*ksound_Hole / 2 
				    / massHole  * (kmag / ksound_Hole 
				    * cos(theta_phonon) - 1 ));

    G4ThreeVector momentum = G4ThreeVector(aTrack.GetMomentumDirection());
    G4ThreeVector orth = momentum.orthogonal();
    G4ThreeVector newDir = momentum.rotate(orth, theta_charge);
    newDir.rotate(aTrack.GetMomentumDirection(), G4UniformRand()*2*pi);

    aParticleChange.ProposeMomentumDirection(newDir);
    aParticleChange.ProposeEnergy(aTrack.GetKineticEnergy()-qEnergy);  
    ResetNumberOfInteractionLengthLeft();    

    return &aParticleChange;

}

G4double LukeScatteringProcess::MakeTheta(G4double& k, G4double& ks){
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
	pValue=G4UniformRand();                                 // *pValMax;//  *(1+2*ks/k+(ks/k)*(ks/k));   
	//need to multiply by unit 's' to make it dimensionless
	pDensity = (cos(theta)-(ks/k))*(cos(theta)-(ks/k))*sin(theta);
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

G4double LukeScatteringProcess::MakePhi(G4double& k,G4double& ks, G4double& theta){


    G4double phi=acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));

    return phi;
  }

G4bool LukeScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
    return((&aPD==DriftingHole::Definition()));
}
