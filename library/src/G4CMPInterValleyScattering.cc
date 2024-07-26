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
// 20170802  Remove MFP calculation; use scattering-rate model
// 20170809  Replace Edelweiss rate with physical (matrix element) model
// 20170821  Use configuration flag to choose Edelweiss vs. physical rate
// 20180831  Change G4CMPInterValleyScattering to use Lin. and Quad. models
// 20190704  Add selection of rate model by name, and material specific
// 20190904  C. Stanford -- Add 50% momentum flip (see G4CMP-168)
// 20190906  Push selected rate model back to G4CMPTimeStepper for consistency
// 20231122  Remove 50% momentum flip (see G4CMP-375)

#include "G4CMPInterValleyScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPInterValleyRate.hh"
#include "G4CMPIVRateQuadratic.hh"
#include "G4CMPIVRateLinear.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include <math.h>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4PhononPolarization.hh"
// Construcor and destructor

G4CMPInterValleyScattering::G4CMPInterValleyScattering()
  : G4CMPVDriftProcess("G4CMPInterValleyScattering", fInterValleyScattering) {
  UseRateModel(G4CMPConfigManager::GetIVRateModel());
}

G4CMPInterValleyScattering::~G4CMPInterValleyScattering() {;}


// Only electrons have physical valleys associated with them

G4bool 
G4CMPInterValleyScattering::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsElectron(&aPD);
}


// Select different rate models by string (globally or by material)

void G4CMPInterValleyScattering::UseRateModel(G4String model) {
  if (model.empty()) {			// No argument, initialize w/global
    if (GetRateModel()) return;		// Do not change existing model
    model = G4CMPConfigManager::GetIVRateModel();
  }

  model.toLower();
  if (model == modelName) return;	// Requested model already in use

  // Select from valid names; fall back to Quadratic if invalid name specified
       if (model(0) == 'q') UseRateModel(new G4CMPIVRateQuadratic);
  else if (model(0) == 'l') UseRateModel(new G4CMPIVRateLinear);
  else if (model(0) == 'i') UseRateModel(new G4CMPInterValleyRate);
  else {
    G4cerr << GetProcessName() << " ERROR: Unrecognized rate model '"
	   << model << "'" << G4endl;
    if (!GetRateModel()) UseRateModel("Quadratic");
  }

  modelName = model;

  // Ensure that TimeStepper process is given new model
  PushModelToTimeStepper();
}


// Switch rate models if necessary based on material

G4double G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& track,
						     G4double prevStep,
						     G4ForceCondition* cond) {
  UseRateModel(theLattice->GetIVModel());	// Use current material's rate
  return G4CMPVProcess::GetMeanFreePath(track, prevStep, cond);
}


// Perform scattering action

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack, 
					 const G4Step& aStep) {
  aParticleChange.Initialize(aTrack); 
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  
  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepDoIt: Step limited by "
	   << postStepPoint->GetProcessDefinedStep()->GetProcessName()
	   << G4endl;
  }
  
  // Don't do anything at a volume boundary
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
//   // Get track's energy in current valley
//   G4ThreeVector p = GetLocalMomentum(aTrack);
//   G4int valley = GetValleyIndex(aTrack);
//   p = theLattice->MapPtoK_valley(valley, p); // p is actually k now
  
//   // picking a new valley at random if IV-scattering process was triggered
//   valley = ChangeValley(valley);
//   G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);

//   p = theLattice->MapK_valleyToP(valley, p); // p is p again
//   RotateToGlobalDirection(p);
    
  G4double density = theLattice->GetDensity();
  G4int nValley = theLattice->NumberOfValleys()-1;		// From symmetry
  G4double m_DOS = theLattice->GetElectronDOSMass();
  G4double m_DOS3half = sqrt(m_DOS*m_DOS*m_DOS);
   G4double alpha = theLattice->GetAlpha();
    G4double Etrk = GetKineticEnergy(aTrack);
    
  std::vector<G4double> probabilities;
    
  G4double scale = nValley*/*kT**/m_DOS3half / (sqrt(2)*pi*hbar_Planck*hbar_Planck*density);

  G4double total = 0.;
  G4int N_op = theLattice->GetNIVDeform();
  for (G4int i = 0; i<N_op; i++) {
    G4double Emin_op = theLattice->GetIVEnergy(i);
    if (Etrk <= Emin_op) {
        probabilities.push_back(0.); 
        continue;		// Apply threshold behaviour
    }

    G4double D_op = theLattice->GetIVDeform(i);
    G4double oscale = scale * D_op*D_op / Emin_op;

    G4double Efunc = sqrt((Etrk-Emin_op)*(1+alpha*(Etrk-Emin_op)))*(1+2*alpha*(Etrk-Emin_op));
      
    G4double orate = oscale * Efunc;

//             G4cout 
// 	     << "phonon rate [" << i << "] " << orate/hertz << " Hz"
// 	     << G4endl;
      
    probabilities.push_back(orate); 
    total += orate;     
  }
    
 std::vector<G4double> cumulatives(probabilities.size());
 for (auto& element : probabilities) {
     element /= total;
     }
 
  cumulatives[0]=probabilities[0];
  G4double test = G4UniformRand();
  G4int Ephononi;
    
      for (G4int i = 1; i<probabilities.size(); i++)  {
cumulatives[i] = probabilities[i-1]+probabilities[i];
      }
    
      for (G4int i = 0; i<cumulatives.size(); i++)  {
          G4cout << cumulatives[i] << "  " << test <<G4endl;
if (test < cumulatives[i])  { Ephononi = i; break;
      }
      }
    
          
          
  // Final state kinematics, generated in accept/reject loop below
  G4double theta_phonon=0, phi_phonon=0, q=0, Ephonon=0, Erecoil=0, costheta=0;
  G4ThreeVector qvec, k_recoil, precoil;	// Outgoing wave vectors
    
  G4int ivalley = GetValleyIndex(aTrack);
    
  G4ThreeVector ktrk(0.);
  ktrk = theLattice->MapV_elToK_HV(ivalley, GetLocalVelocityVector(aTrack));
  G4ThreeVector kdir = ktrk.unit();
  G4double kmag = ktrk.mag();

  Etrk = GetKineticEnergy(aTrack);
  Ephonon = theLattice->GetIVEnergy(Ephononi);
  Erecoil = Etrk - Ephonon;
    
    
         

  costheta=G4UniformRand()*(1-sqrt(Ephonon/Etrk))+sqrt(Ephonon/Etrk);
  phi_phonon=G4UniformRand()*twopi;
    
  q=kmag*costheta-kmag*sqrt(costheta*costheta-Ephonon/Etrk);
    
    

    
  qvec = q*kdir;
  qvec.rotate(kdir.orthogonal(), acos(costheta));
  qvec.rotate(kdir, phi_phonon);

  k_recoil=ktrk-qvec;
  precoil = theLattice->MapK_HVtoP(ivalley, k_recoil);
  RotateToGlobalDirection(precoil);	// Update track in world coordinates
    
    G4cout << "costheta : " << costheta << " Ephononi : " << Ephononi << G4endl;
    
  G4cout << "qvec : " << qvec << " q_mag : " << q << G4endl
      << " Etrk : " << Etrk/eV << " Ephonon : " << Ephonon/eV << " Erecoil : " << Erecoil/eV << G4endl
      << "k_recoil " << k_recoil << " precoil : " << precoil << " ktrk : " << ktrk << G4endl;


    
  // Create real phonon to be propagated, with random polarization
  // If phonon is not created, register the energy as deposited
  G4double weight =
    G4CMP::ChoosePhononWeight(G4CMPConfigManager::GetLukeSampling());
  if (weight > 0.) {
    G4Track* phonon = G4CMP::CreatePhonon(aTrack,
					  G4PhononPolarization::UNKNOWN,
                                          qvec, Ephonon,
                                          aTrack.GetGlobalTime(),
                                          aTrack.GetPosition());
    // Secondary's weight has to be multiplicative with its parent's
    phonon->SetWeight(aTrack.GetWeight() * weight);
    if (verboseLevel>1) {
      G4cout << "phonon wt " << phonon->GetWeight()
	     << " : track " << aTrack.GetTrackID()
	     << " wt " << aTrack.GetWeight()
	     << "  thrown wt " << weight << G4endl;
    }

    aParticleChange.SetSecondaryWeightByProcess(true);
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(phonon);

    // If user wants to track phonons immediately, put track back on stack
    if (secondariesFirst && aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  } else {
    aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);
  }




  // Adjust track kinematics for new valley
  FillParticleChange(ivalley, precoil);
  
  ClearNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}



// Ensure the same rate model is used here and in G4CMPTimeStepper

void G4CMPInterValleyScattering::PushModelToTimeStepper() {
  if (verboseLevel>1)
    G4cout << GetProcessName() << "::PushModelToTimeStepper" << G4endl;

  G4CMPTimeStepper* tsProc =
    dynamic_cast<G4CMPTimeStepper*>(G4CMP::FindProcess(GetCurrentTrack(),
						       "G4CMPTimeStepper"));

  if (tsProc) tsProc->UseIVRateModel(GetRateModel());
}
