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
  
  //If there is no field, there is no IV scattering... but then there
  //is no e-h transport either...
  if (!fMan || !fMan->DoesFieldExist()) return DBL_MAX;

  G4double velocity = GetVelocity(aTrack);
  
  /*  G4double posVec[4] = { 4*0. };
  GetLocalPosition(aTrack, posVec);

  const G4Field* field = fMan->GetDetectorField();
  G4double fieldValue[6];
  field->GetFieldValue(posVec,fieldValue);

  G4ThreeVector fieldVector(fieldValue[3], fieldValue[4], fieldValue[5]);

  if (verboseLevel > 1) {
    G4cout << "IV local position (" << posVec[0] << "," << posVec[1] << ","
	   << posVec[2] << ")\n field " << fieldVector/volt*cm << " V/cm"
	   << "\n magnitude " << fieldVector.mag()/volt*cm << " V/cm toward "
	   << fieldVector.cosTheta() << " z" << G4endl;
	   }

  // Find E-field in HV space: in lattice frame, rotate into valley,
  // then apply HV tansform.
  // NOTE:  Separate steps to avoid matrix-matrix multiplications
  theLattice->RotateToLattice(fieldVector);
  fieldVector *= GetValley(aTrack);
  fieldVector *= theLattice->GetSqrtInvTensor();
  fieldVector /= volt/m;			// Strip units for MFP below

  if (verboseLevel > 1) {
    G4cout << " in HV space " << fieldVector*0.01 << " ("
	   << fieldVector.mag()*0.01 << ") V/cm" << G4endl;
	   }*/

  // Compute mean free path per Edelweiss LTD-14 paper
  // G4double E_0 = theLattice->GetIVField() / (volt/m);
  // G4double mfp = velocity / ( theLattice->GetIVRate() *
  // pow((E_0*E_0 + fieldVector.mag2()), theLattice->GetIVExponent()/2.0) );
  //
  //Acoustic Phonon Scattering 
  G4double energy = GetEnergy(aTrack);
  const  G4double pi = 3.14159265359;
  G4double T =0 ;
  G4double K_b = 1.38064852 *  pow(10,-23);
  G4double h = 1.0545718 * pow (10 , -34);
  G4double M_D = 0;
  G4double D_ac = 0;
  G4double rho = 0; 
  G4double u = 0;
  G4double alpha = 0 
  G4double amfp =(sprt(2)*K_b * T * pow(M_D , 3/2)* (D_ac*D_ac)*
		 (sprt(energy + alpha(enery*energy)))* (1+ (2*alpha*energy)))/
                 (pi* pow(h,4) * rho * ( u*u));
  cout << "this is Acoustic phonon Scattering " <<  amfp << G4endl;

  // Optical phonon Scattering equation 
  G4double D_op = 0 ;
  G4double w_op = 0 ;
  G4double omfp = ( K_b * T * pow (M_D ,3/2) * (D_op*D_op)*
		    (sprt((energy - h*w_op)(1 + alpa(energy - h*w_op))) *
		     (1 + 2*alpha(energy - h*w_op)))) / (sprt(2) * pi*row(h,2)*rho*h*w_op);
  cout << " this is the Optical Phonon Scattering " << omfp << endl ; 
  //Neatral Impurities 
  G4double E_T = 5 * row (10 , -4); 
  G4double m _ = 0 ;
  G4double n_l = 0 ;
   G4double Gamma = (4*sprt(2)* n_l * (h*h) * row(energy,1/2))/
                     (row(m , 3/2)* (energy + E_T));

 cout << " this Neutral Impurities " << Gamma << endl ; 

 G4double mfp  =  velocity / ( Gamma + omfp + amfp ); 

 cout << " this is the mean free path" << mfp << endl ; 
     

  if (verboseLevel > 1) G4cout << "IV MFP = " << mfp/m << G4endl;
  return mfp;
}/*

G4VParticleChange* 
G4CMPIVScatteringPhysical::PostStepDoIt(const G4Track& aTrack, 
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
  }*/

