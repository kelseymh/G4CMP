
#include "PhononScatteringProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "LPhonon.hh"
#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "PhononTrackInformation.hh"
#include "PhysicalLattice.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"


PhononScatteringProcess::PhononScatteringProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{

   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

PhononScatteringProcess::~PhononScatteringProcess()
{;}

PhononScatteringProcess::PhononScatteringProcess(PhononScatteringProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  PhononScatteringProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{

  PhysicalLattice* Lattice = LatticeManager2::getPhysicalLattice(aTrack.GetVolume());
  if(Lattice==0) G4cout<<"\n\nPhononScatteringProcess::PostStepDoIt: WARNING!! PHYSICAL LATTICE POINTER IS NULL!!!\n\n";

  G4double B=Lattice->getScatteringConstant();
  G4double h=6.626e-34*m2*kg/s;
  G4double E= aTrack.GetKineticEnergy();
  
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*B)*aTrack.GetVelocity();


   *condition = NotForced;
   //return DBL_MAX;
   return mfp;
}

G4VParticleChange*
  PhononScatteringProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step& aStep)
{

  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if(postStepPoint->GetStepStatus()==fGeomBoundary)
   { return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);}

  //G4Navigator* theNavigator =
  //    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    aParticleChange.Initialize(aTrack);
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(1);

    // G4cout<<"\nPhononScatteringProcess: Scattering occurred...";

    //randomly generate a new direction
    //modeMixer determines what the new 
    //polarization type will be
    G4Track* sec;
    G4ThreeVector vgroup;  
    G4ThreeVector newDir = G4RandomDirection();
    G4double modeMixer = G4UniformRand();
    
     PhysicalLattice* Lattice = LatticeManager2::getPhysicalLattice(aTrack.GetVolume());
     double cProbST=Lattice->getSTDOS();
     double cProbFT=Lattice->getFTDOS()+cProbST;

    //Generate the new track after scattering
    //the probabilities for the different po-
    //larization types depends on the DOS
    if(modeMixer<cProbST){  
      vgroup=Lattice->mapKtoVDir(1, newDir);
      vgroup=Lattice->LocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(TPhononSlow::PhononDefinition(),vgroup, aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    }else if(modeMixer<cProbFT){
      vgroup=Lattice->mapKtoVDir(2, newDir);
      vgroup=Lattice->LocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(TPhononFast::PhononDefinition(),vgroup, aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    } else {
      vgroup=Lattice->mapKtoVDir(0, newDir);
      vgroup=Lattice->LocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(LPhonon::PhononDefinition(),vgroup, aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    }

    sec->SetUserInformation(new PhononTrackInformation(Lattice->LocalToGlobal.TransformAxis(newDir)));
    aParticleChange.AddSecondary(sec);

    return &aParticleChange;
}

G4bool PhononScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  return ((&aPD==LPhonon::PhononDefinition())|(&aPD==TPhononFast::PhononDefinition())|(&aPD==TPhononSlow::PhononDefinition()));

}
