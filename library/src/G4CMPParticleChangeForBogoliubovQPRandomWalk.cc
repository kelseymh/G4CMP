/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPParticleChangeForBogoliubovQPRandomWalk.cc
/// \brief Implementation of the G4CMPParticleChangeForBogoliubovQPRandomWalk class

#include "G4CMPParticleChangeForBogoliubovQPRandomWalk.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4CMPParticleChangeForBogoliubovQPRandomWalk::G4CMPParticleChangeForBogoliubovQPRandomWalk()
  : G4VParticleChange(){
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4CMPParticleChangeForBogoliubovQPRandomWalk::G4CMPParticleChangeForBogoliubovQPRandomWalk() " << G4endl;
  }
#endif
}

G4CMPParticleChangeForBogoliubovQPRandomWalk::~G4CMPParticleChangeForBogoliubovQPRandomWalk()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4CMPParticleChangeForBogoliubovQPRandomWalk::~G4CMPParticleChangeForBogoliubovQPRandomWalk() " << G4endl;
  }
#endif
}

G4CMPParticleChangeForBogoliubovQPRandomWalk::
G4CMPParticleChangeForBogoliubovQPRandomWalk(const G4CMPParticleChangeForBogoliubovQPRandomWalk &right)
  : G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4CMPParticleChangeForBogoliubovQPRandomWalk::  copy constructor is called " << G4endl;
   }
      theMomentumDirection = right.theMomentumDirection;
      thePosition = right.thePosition;
      theVelocity = right.theVelocity;
}

// assignment operator
G4CMPParticleChangeForBogoliubovQPRandomWalk & G4CMPParticleChangeForBogoliubovQPRandomWalk::operator=(
                                   const G4CMPParticleChangeForBogoliubovQPRandomWalk &right)
{
   if (verboseLevel>1) {
    G4cout << "G4CMPParticleChangeForBogoliubovQPRandomWalk:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
      theTrueStepLength = right.theTrueStepLength;
      theMomentumDirection = right.theMomentumDirection;
      thePosition = right.thePosition;
      theTimeChange = right.theTimeChange;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* G4CMPParticleChangeForBogoliubovQPRandomWalk::UpdateStepForAlongStep(G4Step* pStep)
{
  //  Update the G4Step specific attributes
  pStep->SetStepLength(theTrueStepLength);
  theStatusChange = pStep->GetTrack()->GetTrackStatus();

  // Multiple scattering calculates the final state of the particle
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint->SetPosition( thePosition );
  
  // update velocity
  pPostStepPoint->SetVelocity(theVelocity);
    
  return pStep;
}

G4Step* G4CMPParticleChangeForBogoliubovQPRandomWalk::UpdateStepForPostStep(G4Step* pStep)
{
  // Multiple scattering calculates the final state of the particle
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint->SetPosition( thePosition );
    
  // update velocity
  pPostStepPoint->SetVelocity(theVelocity);
    
  return pStep;
}

//----------------------------------------------------------------
// methods for printing messages
//

void G4CMPParticleChangeForBogoliubovQPRandomWalk::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Position - x (mm)   : "
       << std::setw(20) << thePosition.x()/mm
       << G4endl;
  G4cout << "        Position - y (mm)   : "
       << std::setw(20) << thePosition.y()/mm
       << G4endl;
  G4cout << "        Position - z (mm)   : "
       << std::setw(20) << thePosition.z()/mm
       << G4endl;
  G4cout << "        Momentum Direct - x : "
       << std::setw(20) << theMomentumDirection.x()
       << G4endl;
  G4cout << "        Momentum Direct - y : "
       << std::setw(20) << theMomentumDirection.y()
       << G4endl;
  G4cout << "        Momentum Direct - z : "
       << std::setw(20) << theMomentumDirection.z()
       << G4endl;
  G4cout << "        Velocity : "
         << std::setw(20) << theVelocity
         << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4CMPParticleChangeForBogoliubovQPRandomWalk::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // check

  // MomentumDirection should be unit vector
  accuracy = std::fabs(theMomentumDirection.mag2()-1.0);
  if (accuracy > accuracyForWarning) {
    itsOK = false;
    exitWithError = (accuracy > accuracyForException);
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChangeForMSC::CheckIt  : ";
    G4cout << "the Momentum Change is not unit vector !!"
	   << "  Difference:  " << accuracy << G4endl;
    G4cout << aTrack.GetDefinition()->GetParticleName()
	   << " E=" << aTrack.GetKineticEnergy()/MeV
	   << " pos=" << aTrack.GetPosition().x()/m
	   << ", " << aTrack.GetPosition().y()/m
	   << ", " << aTrack.GetPosition().z()/m
	   <<G4endl;
#endif
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) DumpInfo();
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4CMPParticleChangeForBogoliubovQPRandomWalk::CheckIt",
		"300",
		EventMustBeAborted,
		"momentum direction was illegal");
  }
  //correction
  if (!itsOK) {
    G4double vmag = theMomentumDirection.mag();
    theMomentumDirection = (1./vmag)*theMomentumDirection;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}
