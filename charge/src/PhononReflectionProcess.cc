
#include "PhononReflectionProcess.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "LPhonon.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "PhononTrackInformation.hh"
#include "LatticeManager2.hh"
#include "QPDownconversion.hh"


PhononReflectionProcess::PhononReflectionProcess(const G4String& aName)
:G4VDiscreteProcess(aName)
{
   alminum = NULL;
   kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
   G4cout<<"\n PhononReflectionProcess::Constructor: Geometry surface tolerance is: " << kCarTolerance /mm << " mm";
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

PhononReflectionProcess::~PhononReflectionProcess()
{;}

PhononReflectionProcess::PhononReflectionProcess(PhononReflectionProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  PhononReflectionProcess::GetMeanFreePath( 
       const G4Track& /*aTrack*/, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{
// Always return DBL_MAX and Forced
   *condition = Forced;

   return DBL_MAX;
}


G4VParticleChange*
  PhononReflectionProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step& aStep )
{

  
  PhononTrackInformation* info = (PhononTrackInformation*) (aTrack.GetUserInformation());
  G4ThreeVector reflectedDirection;
  //info->setK(aTrack.GetMomentumDirection());
  

   aParticleChange.Initialize(aTrack);
   G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
   
   // do nothing but return if the current step is not limited by a volume boundary
   if(postStepPoint->GetStepStatus()!=fGeomBoundary)
   { return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep); }
   
     
   //do nothing but return is the step is too short
            if(aTrack.GetStepLength()<=kCarTolerance/2)
   { return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep); }

	    
   const G4DynamicParticle* theDP = aTrack.GetDynamicParticle();
   G4ThreeVector incidentDirection = theDP->GetMomentumDirection();
   G4Navigator* theNavigator = 
     G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
   G4bool valid = true;
   G4ThreeVector localNormal = theNavigator->GetLocalExitNormal(&valid);
   if(valid)
   { localNormal = -localNormal; }
   G4ThreeVector globalNormal = 
     theNavigator->GetLocalToGlobalTransform().TransformAxis(localNormal);

   //Set specular scattering probability SSP//
   G4double SSP = 0;
      
    
   /*
   if(G4UniformRand()<SSP){


     if( incidentDirection*globalNormal>0.0)
       {
	 // this should not happen but .......
	 globalNormal = - globalNormal;;
       }
     G4double PdotN = incidentDirection*globalNormal;

     reflectedDirection = incidentDirection - (2.*PdotN)*globalNormal;
     info->setK(reflectedDirection);  //This is just a temporary bug fix, in order to determine a k-vector. Not physical
     
     //reflectedDirection=LatticeManager2::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);

   } else {
     /////////If scattering is diffuse:///////////
     G4ThreeVector reflectedK;

     G4double PdotN = incidentDirection*globalNormal;
     if(PdotN>0.)
       {
	 // this should not happen but .......
	 globalNormal = - globalNormal;
	 PdotN *= -1.;
       }
     //find reflected direction according to lambert's rule
     //G4ThreeVector aux=G4RandomDirection();
     //while(globalNormal.cross(aux).mag()==0) aux=G4RandomDirection();
     //reflectedDirection = (globalNormal.rotate(globalNormal.cross(aux), acos(2*G4UniformRand()))*0.90);
     reflectedDirection = G4LambertianRand(globalNormal);
   
     info->setK(reflectedDirection); 
     //if(aTrack.GetDefinition()==LPhonon::PhononDefinition()) reflectedDirection=LatticeManager2::mapKtoVDir(aTrack.GetVolume(),0,reflectedDirection);
     //else if(aTrack.GetDefinition()==TPhononSlow::PhononDefinition()) reflectedDirection=LatticeManager2::mapKtoVDir(aTrack.GetVolume(),1,reflectedDirection);
     //else if(aTrack.GetDefinition()==TPhononFast::PhononDefinition()) reflectedDirection=LatticeManager2::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);

   }



   aParticleChange.ProposeMomentumDirection(reflectedDirection.unit());

   //   check if phonon is lost to black body radiation
   if(postStepPoint->GetMaterial()!=alminum)
     if(G4UniformRand()<0.001) { 
       aParticleChange.ProposeTrackStatus(fStopAndKill);
   }
   // in case the other side is alminum electrode, deposit energy
   //   QPDownconversion con;

        if((postStepPoint->GetMaterial()==alminum)&&(G4UniformRand()<0.02013)&&(aTrack.GetKineticEnergy()>347.43e-6*eV))
     //if(!LatticeManager2::hasLattice(postStepPoint->GetPhysicalVolume()))
     {
       //       con.downconvert(aTrack, &aParticleChange, reflectedDirection);
       G4double TwoAlGap = 347.43e-6*eV;
       G4double eKin = aTrack.GetKineticEnergy();
     
       if(eKin<4*TwoAlGap){
	 aParticleChange.ProposeTrackStatus(fStopAndKill);   
	  aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
       } else{
        
	 aParticleChange.ProposeNonIonizingEnergyDeposit(4*TwoAlGap);
	 aParticleChange.ProposeEnergy(eKin-4*TwoAlGap);
	 //if(G4UniformRand()>(5.6/15.0)){
	 //aParticleChange.ProposeTrackStatus(fStopAndKill);
	 //aParticleChange.ProposeEnergy(0);
	   //}
      
       }
       
     }else if(aTrack.GetKineticEnergy()<347.4e-6*eV){

	 aParticleChange.ProposeTrackStatus(fStopAndKill);  
     }
   */
	
     //This stops ANY phonon encountering a wall. Debugging purposes.
	  if(postStepPoint->GetMaterial()==alminum){
	  G4double eKin = aTrack.GetKineticEnergy();
	  aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
	  aParticleChange.ProposeTrackStatus(fStopAndKill);
	  }

	  aParticleChange.ProposeTrackStatus(fStopAndKill);
	
   return &aParticleChange;
}

G4bool PhononReflectionProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  return ((&aPD==TPhononFast::PhononDefinition())||(&aPD==LPhonon::PhononDefinition())||(&aPD==TPhononSlow::PhononDefinition()));
}

void PhononReflectionProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
   if(!alminum)
   { alminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"); }
}

