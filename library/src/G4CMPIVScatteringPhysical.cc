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
git}



G4double 
G4CMPIVScatteringPhysical::GetMeanFreePath(const G4Track& aTrack,
					    G4double,
					    G4ForceCondition* condition) {
 *condition = NotForced;
 

  G4double velocity = GetVelocity(aTrack);
  G4double energy = GetKineticEnergy(aTrack);
  //Acoustic Phonon Scattering 

  G4double T = 0.015*kelvin ;
  G4double mass_electron = electron_mass_c2/c_squared;  
       // 9.109e-31*kilogram ;	// killograms?  electron mass? same as M_D
 
  G4double ml = 1.38 *mass_electron;
  G4double mt = .081 *mass_electron;

  G4double M_D =cbrt( ml*mt*mt)* kilogram ;	// kg  density of mass state get from Lattice class (eventually) 
  G4double D_ac = 11*eV;            // 1.7622e-18;/ Units?  defermation material
  G4double rho = 5.327e3 * kilogram/m3 ; 		
  G4double alpha=0.3/eV;                         // 1.872659176e18  * 1/coulomb ;     

  //Units?  freespace?
  G4double epsilon_r = 16.2;
  G4double epsilon = epsilon_r * epsilon0;	      

  // Useful constants for expressions below
  const G4double hbar_sq  = hbar_Planck * hbar_Planck;
  const G4double hbar_4th = hbar_sq * hbar_sq;
  const G4double M_D3half = sqrt(M_D*M_D*M_D);
  const G4double D_ac_sq = D_ac*D_ac;
  const G4double alpha_times_energy = alpha*energy;
  const G4double velocity_sq = velocity*velocity;

  G4double amfp =(sqrt(2)*k_Boltzmann * T * M_D3half * D_ac_sq *
		  sqrt(energy + alpha_times_energy *energy) * (1+ 2*alpha_times_energy))/
                 (pi* hbar_4th * rho * velocity_sq);
  

  // Optical phonon Scattering equation
  G4double D_op []= { 3e10*eV, 2e9*eV }  ;
  G4double w_op []= { 27.3e-3*eV, 10.3e-3*eV } ;
  G4double omfp[] = {0,0};
  G4double omfpTotal = 0;
  
  for (int i = 0; i<2; i++) {
    G4double hw_op = hbar_Planck * w_op[i];// Energy of optical Phonon
    G4double energy_minus_hw_op = energy - hw_op;
    G4double D_op_sq = D_op[i]*D_op[i];
    G4double alpha_times_ehw_op = alpha * energy_minus_hw_op;
    G4double everything_under_sqrt = sqrt(energy_minus_hw_op + energy_minus_hw_op * alpha_times_ehw_op);

    omfp[i] =  k_Boltzmann * T * M_D3half * D_op_sq * everything_under_sqrt *
      (1 + 2*alpha_times_ehw_op) / (sqrt(2) * pi* hbar_sq *rho*hw_op);
   
    omfpTotal+=omfp[i];
  }
 
  //Neutral Impurities
  G4double e_r =  epsilon0/epsilon;
  G4double E_T = (M_D /mass_electron) * e_r;
  G4double n_l = 1e17 ;			// Units? The number density of inpurities
  
  G4double Gamma = (4*sqrt(2)* n_l * hbar_sq * sqrt(energy))/
    ( sqrt(m*m*m)* (energy + E_T));
  
  
  
  G4double mfp  =  velocity / ( Gamma + omfpTotal + amfp ); 
  
  
  
  if (verboseLevel > 1) {
    G4cout << "IV MFP = " << mfp/m << G4endl;
    G4cout << "this is Acoustic phonon Scattering " <<  amfp << G4endl;
    // G4cout << " this is the Optical Phonon Scattering " << i << " " << ofmp[i] << G4endl;   
    G4cout << " this Neutral Impurities " << Gamma << G4endl ;
    G4cout << " this is the mean free path" << mfp << G4endl ;       
    return mfp; }
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

