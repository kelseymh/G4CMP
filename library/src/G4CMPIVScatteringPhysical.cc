/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140324  Drop hard-coded IV scattering parameters; get from lattice
// 20140324  Restore Z-axis mass tensor
// 20140331  Add required process subtype code
// 20140418  Drop local valley transforms, use lattice functions instead
// 20140429  Recompute kinematics relative to new valley
// 20140908  Allow IV scatter to change momentum by conserving energy
// 20150109  Revert IV scattering to preserve momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160601  Must apply lattice rotation before valley.
// 20161004  Change valley selection function to avoid null choice
// 20161114  Use G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for particle identity checks

//#include "G4CMPInterValleyScattering.hh"
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
#include "G4CMPIVScatteringPhysical.hh"
G4CMPIVScatteringPhysical::G4CMPIVScatteringPhysical()
  : G4CMPVDriftProcess("G4CMPIVScatteringPhysical", fIVScatteringPhysical) {;}

G4CMPIVScatteringPhysical::~G4CMPIVScatteringPhysical() {;}


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
  const  G4double pi = 3.14159265359;
  G4double T = 0.015 ;
  G4double K_b = 1.38064852e-23;
  G4double h = 1.0545718e-34;
  G4double M_D = 1.98615472 * pow (10, -31) ;
  G4double D_ac = 1.7622* pow( 10, -18);
  G4double rho = 5.327 * pow(10 , 3); 
  G4double alpha = 1.872659176 * pow(10 ,18);
  G4double m = 9.109 * pow(10, -31);
  G4double e_0 = 8.85 * pow(10, -12);
  G4double e = 1.4337 * pow(10, -10);

  G4double amfp =(sqrt(2)*K_b * T * pow(M_D , 1.5)* (D_ac*D_ac)*
		 (sqrt(energy + alpha(enery*energy)))* (1+ (2*alpha*energy)))/
                 (pi* pow(h,4) * rho * ( velocity*velocity));
  cout << "this is Acoustic phonon Scattering " <<  amfp << G4endl;

  // Optical phonon Scattering equation 
  G4double D_op []= {3 * pow(10,10),2* pow(10,9)}  ;
  G4double w_op []= {27.3,10.3} ;
  G4double amfp =(sqrt(2)*K_b * T * pow(M_D , 1.5)* (D_ac*D_ac)*
		  (sqrt(energy + alpha(enery*energy)))* (1+ (2*alpha*energy)))/
    (pi* pow(h,4) * rho * ( velocity*velocity));
  cout << "this is Acoustic phonon Scattering " <<  amfp << G4endl;

  // Optical phonon Scattering equation
  G4double D_op []= {3 * pow(10,10),2* pow(10,9)}  ;
  G4double w_op []= {27.3,10.3} ;
  G4double omfp[] = {0,0};
  G4double omfpTotal = 0;
  for (int i = 0 ;i<2 ; i++){
    D_op[i]=D_op[i]*ev;
    w_op[i]= w_op* .001 * ev;
    omfp[i]=( K_b * T * pow (M_D ,1.5) * (D_op[i]*D_op[i]))*
      sqrt((energy -( h*w_op[i]))*(1 + alpha*(energy -( h*w_op[i])))) *
      (1 + 2*alpha*(energy - h*w_op[i])) / (sqrt(2) * pi* pow(h,2)	\
					    *rho*h*w_op[i]);
    cout << " this is the Optical Phonon Scattering " << i << " " << ofmp[i] << endl;
    omfpTotal+=omfp[i];
  }
 
  //Neutral Impurities 
  G4double E_T =(M_D /m) * (e_o /e) ;
  G4double n_l = pow(10,17) ;
  G4double Gamma = (4*sprt(2)* n_l * (h*h) * row(energy,1/2))/
    (row(m , 3/2)* (energy + E_T));
  
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

