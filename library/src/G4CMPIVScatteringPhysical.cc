/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170711  New class implementing physically motivated IV scattering

#include "G4CMPIVScatteringPhysical.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "math.h"


// Constructor and destructor

G4CMPIVScatteringPhysical::G4CMPIVScatteringPhysical()
  : G4CMPVDriftProcess("G4CMPIVScatteringPhysical", fInterValleyScattering) {;}

G4CMPIVScatteringPhysical::~G4CMPIVScatteringPhysical() {;}


// Flag which kind of particle does IV scattering

G4bool 
G4CMPIVScatteringPhysical::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsElectron(&aPD);
}


G4double 
G4CMPIVScatteringPhysical::GetMeanFreePath(const G4Track& aTrack,
					    G4double,
					    G4ForceCondition* condition) {
  *condition = NotForced;

  // Get electric field associated with current volume, if any
  G4FieldManager* fMan =
    aTrack.GetVolume()->GetLogicalVolume()->GetFieldManager();
  
  //is no e-h transport either...
  if (!fMan || !fMan->DoesFieldExist()) 
    return DBL_MAX;

  G4double velocity = GetVelocity(aTrack);

  //Acoustic Phonon Scattering 
  G4double energy = GetEnergy(aTrack);
  G4double T = 0.015*kelvin ;

  G4double M_D = 1.98615472e-31;	// Units?  Meaning?
  G4double D_ac = 1.7622e-18;		// Units?  Meaning?
  G4double rho = 5.327e3 ; 		// Units?  Meaning?
  G4double alpha = 1.872659176e18 ;	// Units?  Meaning?
  G4double m = 9.109e-31;		// Units?  Meaning?
  G4double e_0 = 8.85e-12 ;		// Units?  Meaning?
  G4double e = 1.4337e-10 ;		// Units?  Meaning?

  // Useful constants for expressions below
  const G4double h_sq  = h_Planck*h_Planck;
  const G4double h_4th = h_sq * h_sq;
  const G4double M_D3half = sqrt(M_D*M_D*M_D);

  G4double amfp =(sqrt(2)*k_Boltzmann * T * M_D3half * D_ac*D_ac *
		 (sqrt(energy + alpha*(energy*energy)))* (1+ (2*alpha*energy)))/
                 (pi* h_4th * rho * velocity*velocity);
  cout << "this is Acoustic phonon Scattering " <<  amfp << G4endl;

  // Optical phonon Scattering equation
  G4double D_op []= { 3e10*eV, 2e9*eV }  ;
  G4double w_op []= { 27.3e-3*eV, 10.3e-3*eV } ;
  G4double omfp[] = {0,0};
  G4double omfpTotal = 0;
  for (int i = 0; i<2; i++) {
    G4double hw_op = h_Planck*w_op[i];

    omfp[i] = ( k_Boltzmann * T * M_D3half * (D_op[i]*D_op[i]))*
      sqrt((energy - hw_op)*(1 + alpha*(energy - hw_op))) *
      (1 + 2*alpha*(energy - hw_op)) / (sqrt(2) * pi* h_sq *rho*hw_op);
    cout << " this is the Optical Phonon Scattering " << i << " " << ofmp[i] << endl;
    omfpTotal+=omfp[i];
  }
 
  //Neutral Impurities 
  G4double E_T = (M_D /m) * (e_0 /e) ;
  G4double n_l = 1e17 ;			// Units?
  G4double Gamma = (4*sqrt(2)* n_l * (h_Planck*h_Planck) * sqrt(energy))/
    ( sqrt(m*m*m)* (energy + E_T));
  
  cout << " this Neutral Impurities " << Gamma << endl ; 
  
  G4double mfp  =  velocity / ( Gamma + omfpTotal + amfp ); 
  
  cout << " this is the mean free path" << mfp << endl ;      
  
  if (verboseLevel > 1) 
    G4cout << "IV MFP = " << mfp/m << G4endl;
  
  return mfp;
}


G4VParticleChange*  G4CMPIVScatteringPhysical::PostStepDoIt(const G4Track& aTrack, 
							    const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  // Do nothing when step limit reached or leaving volume
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by process "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }
  
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return &aParticleChange;
  }
  
  // Get track's energy in current valley
  G4ThreeVector p = GetLocalMomentum(aTrack);
  G4int valley = GetValleyIndex(aTrack);
  p = theLattice->MapPtoK_valley(valley, p); // p is actually k now
  
  // picking a new valley at random if IV-scattering process was triggered
  valley = ChangeValley(valley);
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);
  
  p = theLattice->MapK_valleyToP(valley, p); // p is p again
  RotateToGlobalDirection(p);
  
  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);
  
  ResetNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}

