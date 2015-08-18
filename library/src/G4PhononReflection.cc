//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/src/G4PhononReflection.cc
/// \brief Implementation of the G4PhononReflection class
//
// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly 
// geometry dependent.Currently, phonons are killed when
// they reach a boundary. If the other side of the 
// boundary was Al, a hit is registered.
//  
// $Id$
//
// 20131115  Throw exception if track's polarization state is invalid.
// 20140103  Move charge version of code here, still commented out
// 20140331  Add required process subtype code

#include "G4CMPConfigManager.hh"
#include "G4PhononReflection.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VParticleChange.hh"
#include "G4LogicalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "Randomize.hh"


G4PhononReflection::G4PhononReflection(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4PhononReflection::~G4PhononReflection() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Always return DBL_MAX and Forced. This ensures that the process is
// called at the end of every step. In PostStepDoIt the process
// decides whether the step encountered a volume boundary and a
// reflection should be applied

G4double G4PhononReflection::GetMeanFreePath(const G4Track&, G4double,
					     G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly geometry
// dependent.Currently, phonons are killed when they reach a
// boundary. If the other side of the boundary was Al, a hit is
// registered.
  
G4VParticleChange* G4PhononReflection::PostStepDoIt(const G4Track& aTrack,
						    const G4Step& aStep) { 
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  G4StepPoint* preStepPoint = aStep.GetPreStepPoint();

  if (postStepPoint->GetStepStatus()!=fGeomBoundary) {
    //make sure that correct phonon velocity is used after the step
    int pol = GetPolarization(aTrack);
    if (pol < 0 || pol > 2) {
      G4Exception("G4PhononReflection::PostStepDoIt","Phonon001",
      EventMustBeAborted, "Track is not a phonon");
      return &aParticleChange;		// NOTE: Will never get here
    }

    // FIXME:  This should be using wave-vector, shouldn't it?
    G4double vg = theLattice->MapKtoV(pol, trackKmap->GetK(aTrack));

    //Since step was not a volume boundary, just set correct phonon velocity and return
    aParticleChange.ProposeVelocity(vg);
    return &aParticleChange;
  }

  // do nothing but return if the step is too short
  // This is to allow actual reflection where after
  // the first boundary crossing a second, infinitesimal
  // step occurs crossing back into the original volume
  if (aTrack.GetStepLength()<=kCarTolerance/2) {
    return &aParticleChange;
  }

  if (verboseLevel) {
    G4cout << GetProcessName() << "::PostStepDoIt length "
     << aTrack.GetStepLength() << G4endl;
  }

  G4VPhysicalVolume* thePrePV = preStepPoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = postStepPoint->GetPhysicalVolume();

  if (thePrePV == thePostPV) {
    if (verboseLevel) {
      G4cerr << GetProcessName() << " ERROR: fGeomBoundary status set, but"
       << " pre- and post-step volumes are identical!" << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  // do nothing if the current step is inbound from outside the original volume
  G4LatticePhysical* volLattice =
    G4LatticeManager::GetLatticeManager()->GetLattice(thePrePV);
  if (volLattice != theLattice) {
    if (verboseLevel>1) {
      G4cout << GetProcessName() << ": Track inbound after reflection"
       << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  if (verboseLevel>1) {
    G4cout <<   "    Track volume: " << aTrack.GetVolume()->GetName()
     << "\n  PreStep volume: " << thePrePV->GetName() << " @ "
     << preStepPoint->GetPosition()
     << "\n PostStep volume: " << thePostPV->GetName() << " @ "
     << postStepPoint->GetPosition()
     << G4endl;
  }

  // Grab surface information
  G4LogicalSurface* surface = G4LogicalBorderSurface::GetSurface(thePrePV,
                                                                 thePostPV);
  G4CMPSurfaceProperty* borderSurface;
  if (surface) {
    borderSurface =
      static_cast<G4CMPSurfaceProperty*>(surface->GetSurfaceProperty());
  } else {
    if (verboseLevel>1) {
      G4cerr << GetProcessName() << ": No border surface defined for "
       << thePrePV->GetName() << " to "  << thePostPV->GetName()
       << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  G4MaterialPropertiesTable* phonPropTable =
    const_cast<G4MaterialPropertiesTable*>(
      borderSurface->GetPhononMaterialPropertiesTablePointer());
  absProb = phonPropTable->GetConstProperty("absProb");
  specProb = phonPropTable->GetConstProperty("specProb");
  gapEnergy = phonPropTable->GetConstProperty("gapEnergy");

  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  surfNorm =
    iNav[navID]->GetGlobalExitNormal(postStepPoint->GetPosition(), &goodNorm);
  if (!goodNorm) {
    G4cerr << GetProcessName() << " ERROR:  Cannot get normal at surface of "
     << thePrePV->GetName() << " @ " << postStepPoint->GetPosition()
     << G4endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  } else if (verboseLevel>2) {
    G4cout << " Normal " << surfNorm << " @ " << postStepPoint->GetPosition()
     << G4endl;
  }

  const G4int maxRefl = G4CMPConfigManager::GetMaxPhononBounces();
  if (maxRefl<0 || numberOfReflections < maxRefl) {
    if (ReflectTrack(aStep)) {
      numberOfReflections++;
      return DoReflection(aStep);
    } else {
      if (verboseLevel) {
        G4cout << GetProcessName() << " WARNING: Phonon has reflected "
          << maxRefl << " times. Track being killed." << G4endl;
      }
      G4double eKin = aTrack.GetKineticEnergy();
      aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }
  }


  /*** THIS IS WHERE REFLECTION/TRANSMISSION CALCULATIONS BELONG ***

  const G4DynamicParticle* theDP = aTrack.GetDynamicParticle();
  G4ThreeVector incidentDirection = theDP->GetMomentumDirection();
  G4Navigator* theNavigator = 
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  G4bool valid = true;
  G4ThreeVector localNormal = theNavigator->GetLocalExitNormal(&valid);
  if (valid) { localNormal = -localNormal; }
  G4ThreeVector globalNormal = 
    theNavigator->GetLocalToGlobalTransform().TransformAxis(localNormal);
  
  //Set specular scattering probability SSP//
  G4double SSP = 0;

  if (G4UniformRand()<SSP){
    if (incidentDirection*globalNormal>0.0) {
      // this should not happen but .......
      globalNormal = - globalNormal;
    }

    G4double PdotN = incidentDirection*globalNormal;

    reflectedDirection = incidentDirection - (2.*PdotN)*globalNormal;
    info->setK(reflectedDirection);  //This is just a temporary bug fix, in order to determine a k-vector. Not physical
    //reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);
  } else {
    /////////If scattering is diffuse:///////////
    G4ThreeVector reflectedK;
    
    G4double PdotN = incidentDirection*globalNormal;
    if (PdotN>0.) {
      // this should not happen but .......
      globalNormal = - globalNormal;
      PdotN *= -1.;
    }
    
    reflectedDirection = G4LambertianRand(globalNormal);
    
    info->setK(reflectedDirection); 
    //if(aTrack.GetDefinition()==G4PhononLong::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),0,reflectedDirection);
    //else if(aTrack.GetDefinition()==G4PhononTransSlow::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),1,reflectedDirection);
    //else if(aTrack.GetDefinition()==G4PhononTransFast::PhononDefinition()) reflectedDirection=G4LatticeManager::mapKtoVDir(aTrack.GetVolume(),2,reflectedDirection);
  }
  
  aParticleChange.ProposeMomentumDirection(reflectedDirection.unit());
  
  //   check if phonon is lost to black body radiation
  if (postStepPoint->GetMaterial() != alminum) {
    if (G4UniformRand()<0.001) { 
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }
  } else {
    // in case the other side is alminum electrode, deposit energy
    //   QPDownconversion con;
    if ((G4UniformRand()<0.02013) && 
	(aTrack.GetKineticEnergy()>347.43e-6*eV)) {
      //       con.downconvert(aTrack, &aParticleChange, reflectedDirection);
      G4double TwoAlGap = 347.43e-6*eV;
      G4double eKin = aTrack.GetKineticEnergy();
      
      if(eKin<4*TwoAlGap){
	aParticleChange.ProposeTrackStatus(fStopAndKill);   
	aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
      } else {
	aParticleChange.ProposeNonIonizingEnergyDeposit(4*TwoAlGap);
	aParticleChange.ProposeEnergy(eKin-4*TwoAlGap);
	//if(G4UniformRand()>(5.6/15.0)){
	//aParticleChange.ProposeTrackStatus(fStopAndKill);
	//aParticleChange.ProposeEnergy(0);
	//}
      }
    } else if (aTrack.GetKineticEnergy()<347.4e-6*eV){
      aParticleChange.ProposeTrackStatus(fStopAndKill);  
    }
  } else
  ***/

  return &aParticleChange; 
}

G4bool G4PhononReflection::ReflectTrack(const G4Step& aStep)
{
  // Track is outbound hitting surface of volume
  return (aStep.GetPostStepPoint()->GetMomentumDirection()*surfNorm > 0.);
}

G4VParticleChange* G4PhononReflection::DoReflection(const G4Step& aStep)
{
  // If the phonon has energy < 2*GapEnergy, we'll never detect it since it
  // can't break a Cooper pair.
  G4double eKin = aStep.GetPostStepPoint()->GetKineticEnergy();
  if (eKin < 2.0*gapEnergy) {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector k = trackKmap->GetK(aStep.GetTrack());
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << k.unit() << G4endl;

  if (G4UniformRand()<specProb){
    G4ThreeVector momDir = k.unit();
    // Specular reflecton reverses momentum along normal
    G4double momNorm = momDir * surfNorm;
    momDir -= 2.*momNorm*surfNorm;

    if (verboseLevel>2)
      G4cout << " New momentum direction " << momDir << G4endl;

    k = k.mag()*momDir;
    trackKmap->SetK(aStep.GetTrack(),k);
    aParticleChange.ProposeMomentumDirection(
      theLattice->MapKtoVDir(GetPolarization(aStep.GetTrack()),k));
  } else {
    G4ThreeVector reflectedkDir = LambertRotation();

    if (verboseLevel>2)
      G4cout << " New momentum direction " << reflectedkDir << G4endl;

    k = k.mag()*reflectedkDir;
    trackKmap->SetK(aStep.GetTrack(),k);
    aParticleChange.ProposeMomentumDirection(
      theLattice->MapKtoVDir(GetPolarization(aStep.GetTrack()),k));
  }

  //We know eKin >= 2.0*gapEnergy
  aParticleChange.ProposeNonIonizingEnergyDeposit(eKin - 2.0*gapEnergy);
  if (eKin - 2.0*gapEnergy < 2.0*gapEnergy) {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }
  return &aParticleChange;
}

G4ThreeVector G4PhononReflection::LambertRotation()
{
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1) / 2.0;

  G4ThreeVector refl = -surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}
