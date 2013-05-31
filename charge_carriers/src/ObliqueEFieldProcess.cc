 #include "ObliqueEFieldProcess.hh"
 #include "G4Step.hh"
 #include "G4Track.hh"
 #include "G4ThreeVector.hh"
 #include <fstream>
 #include <iostream>

 ObliqueEFieldProcess::ObliqueEFieldProcess(const G4String& aName) : G4VDiscreteProcess(aName)
 {
   normalToValley= G4AffineTransform(G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532));
   valleyToNormal= G4AffineTransform(G4RotationMatrix(G4ThreeVector(1.0,0,0), 0.95532)).Inverse();

  if(verboseLevel>1){
    G4cout<<GetProcessName()<<" is created "<<G4endl;
  }
}

ObliqueEFieldProcess::~ObliqueEFieldProcess()
{ ; }

ObliqueEFieldProcess::ObliqueEFieldProcess(ObliqueEFieldProcess& right)
  : G4VDiscreteProcess(right)
{ ; }

G4double ObliqueEFieldProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  *condition=Forced;
  return DBL_MAX;
}

G4VParticleChange* ObliqueEFieldProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  G4double e = -1.0*eplus;
  G4double me=9.1e-31*kg;
  G4double mx=0.081*me;
  G4double my=0.081*me;
  G4double mz=1.58*me;
  G4double EField= 20.0*volt/m;
  G4ThreeVector EDir = G4ThreeVector(0.0, 0.0, 1.0).unit();

  //calculate oblique momentum change
  EDir = normalToValley.TransformAxis(EDir).unit();
  double dvx = e*EField*aStep.GetDeltaTime()*EDir.getX()/mx;
  double dvy = e*EField*aStep.GetDeltaTime()*EDir.getY()/my;
  double dvz = e*EField*aStep.GetDeltaTime()*EDir.getZ()/mz;
  G4ThreeVector deltaV = G4ThreeVector(dvx, dvy, dvz).unit();
  deltaV.setMag(sqrt(dvx*dvx + dvy*dvy + dvz*dvz));
  deltaV = valleyToNormal.TransformAxis(deltaV);

  G4ThreeVector oldV=aTrack.GetMomentumDirection();
  oldV.setMag(aTrack.GetVelocity());

  G4double dt = aStep.GetDeltaTime();
  G4ThreeVector posn = aStep.GetPostStepPoint()->GetPosition();
  G4ThreeVector dx = G4ThreeVector(deltaV.getX()*dt,
					deltaV.getY()*dt,
					deltaV.getY()*dt);
  G4ThreeVector newPosn = G4ThreeVector(posn.getX() + dx.getX(),
					posn.getY() + dx.getY(),
					posn.getZ() + dx.getZ()
					);

  //calculate new momentum  
  G4ThreeVector newV;
  newV.set(oldV.getX() + deltaV.getX(),oldV.getY() + deltaV.getY(),oldV.getZ() + deltaV.getZ()) ;

  G4double dE=0.5*0.12*9.1e-31*((newV.mag()*newV.mag()) - (oldV.mag()*oldV.mag()))/(m/s)/(m/s)/1.602e-19*eV;
  G4double newE=aTrack.GetKineticEnergy()+dE;

  G4double tim = aStep.GetPostStepPoint()->GetGlobalTime();

  /*
  G4cout<<"\n--------------------------";
  G4cout<<"\nEFieldProcess::PostStepDoIt: dt: " <<aStep.GetDeltaTime()/s;
  G4cout<<"\n\tcorresponding to step length dx="<<(aTrack.GetVelocity()*aStep.GetDeltaTime())/mm <<" mm";
  G4cout<<"\nEFieldProcess::PostTepDoIt:dV: "<<dv/(m/s);
  G4cout<<"\n\tdv.x:"<<deltaV.getX()/(m/s);
  G4cout<<"\n\tdv.y:"<<deltaV.getY()/(m/s);
  G4cout<<"\n\tdv.z:"<<deltaV.getZ()/(m/s);
  G4cout<<"\nEFieldProcess::PostStepDoIt: Previous velocity oldV: "<<aTrack.GetVelocity()/(m/s);
  G4cout<<"\n\toldV(x,y,z): ("<<oldV.getX()/(m/s)<<","<<oldV.getY()/(m/s)<<","<<oldV.getZ()/(m/s)<<")";
  G4cout<<"\nEFieldProcess::PostStepDoIt: new velocity newV: "<<newV.mag()/(m/s);
  G4cout<<"\n\tnewV(x,y,z): ("<<newV.getX()/(m/s)<<","<<newV.getY()/(m/s)<<","<<newV.getZ()/(m/s)<<")";
  G4cout<<"\ncorresponding to acceleration a="<<(newV.mag()-aTrack.GetVelocity())/(aStep.GetDeltaTime())/(m/s/s);
  G4cout<<"\nEFieldProcess::PostStepDoIt: previous kE = "<<aTrack.GetKineticEnergy()/eV;
  G4cout<<"\n\tnew kE: "<< newE/eV;
  G4cout<<"difference in energies: "<<(newE/eV - aTrack.GetKineticEnergy()/eV);
  
  G4cout<<"\n";
  */
  aParticleChange.Initialize(aTrack);

  aParticleChange.ProposeMomentumDirection(newV.unit());
  aParticleChange.ProposeEnergy(newE);
  //aParticleChange.ProposeGlobalTime(tim + 0.1*dt);
  aParticleChange.ProposePosition(newPosn);

  std::ofstream epositions;
  epositions.open("e-positions.txt", std::ofstream::app);

  epositions << newPosn.getX() << " " << newPosn.getY() << " " << newPosn.getZ() << "\n";
  epositions.close();

  return &aParticleChange;
}

G4bool ObliqueEFieldProcess::IsApplicable(const G4ParticleDefinition& aPD){
  return true;
}
