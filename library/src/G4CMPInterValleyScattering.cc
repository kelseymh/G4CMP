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
// 20240823  Allow ConfigManager IVRateModel setting to override config.txt
// 20250424  Add phonon emission and angular distribution.

#include "G4CMPInterValleyScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPInterValleyRate.hh"
#include "G4CMPIVRateQuadratic.hh"
#include "G4CMPIVRateLinear.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"
#include <math.h>


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

    model = (G4CMPConfigManager::GetIVRateModel().empty() ? "Quadratic"
	     : G4CMPConfigManager::GetIVRateModel());
  }

  model.toLower();
  if (model == modelName) return;	// Requested model already in use

  // Select from valid names; fall back to Quadratic if invalid name specified
       if (model(0) == 'q') {
           UseRateModel(new G4CMPIVRateQuadratic);
           doValleySwitch=true;    // Deprecated IV scattering PostStepDoIt
           }
  else if (model(0) == 'l') {
      UseRateModel(new G4CMPIVRateLinear);
      doValleySwitch=true;    // Deprecated IV scattering PostStepDoIt
      }
  else if (model(0) == 'i') {
      UseRateModel(new G4CMPInterValleyRate);
      doValleySwitch=false;    // Up-to-date IV scattering PostStepDoIt
      }
  else {
    G4cerr << GetProcessName() << " ERROR: Unrecognized rate model '"
	   << model << "'" << G4endl;
    if (!GetRateModel()) UseRateModel("Quadratic");
  }

  modelName = GetRateModel()->GetName();
  modelName.toLower();

  // Ensure that TimeStepper process is given new model
  PushModelToTimeStepper();
}


// Switch rate models if necessary based on material

G4double G4CMPInterValleyScattering::GetMeanFreePath(const G4Track& track,
						     G4double prevStep,
						     G4ForceCondition* cond) {
  // If user set a model in ConfigManager, use that
  G4String userModel = G4CMPConfigManager::GetIVRateModel();
  if (!userModel.empty()) UseRateModel(userModel);
  else UseRateModel(theLattice->GetIVModel());	// Use current material's rate

  return G4CMPVProcess::GetMeanFreePath(track, prevStep, cond);
}


// Perform scattering action

G4VParticleChange* 
G4CMPInterValleyScattering::PostStepDoIt(const G4Track& aTrack,
                                       const G4Step& aStep) {
return (doValleySwitch ? SwitchValleys(aTrack, aStep) : ValleyScattering(aTrack, aStep));
}


// Deprecated IV scattering

G4VParticleChange* 
G4CMPInterValleyScattering::SwitchValleys(const G4Track& aTrack, 
					 const G4Step& aStep) {
  InitializeParticleChange(GetValleyIndex(aTrack), aTrack);
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
  
  // Get track's energy in current valley
  G4ThreeVector p = GetLocalMomentum(aTrack);
  G4int valley = GetValleyIndex(aTrack);
  p = theLattice->MapPtoK(valley, p); // p is actually k now
  p = theLattice->RotateToValley(valley, p);
  
  // picking a new valley at random if IV-scattering process was triggered
  valley = ChangeValley(valley);
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);

  p = theLattice->RotateFromValley(valley, p);
  p = theLattice->MapKtoP(valley, p); // p is p again
  RotateToGlobalDirection(p);
  
  // Adjust track kinematics for new valley
  FillParticleChange(valley, p);
  
  ClearNumberOfInteractionLengthLeft();    
  return &aParticleChange;
}


// Up-to-date IV scattering

G4VParticleChange* 
G4CMPInterValleyScattering::ValleyScattering(const G4Track& aTrack,
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

  // Get track's kinematics in current valley
  G4int valley = GetValleyIndex(aTrack);
  G4ThreeVector ptrk = GetLocalMomentum(aTrack);
  G4ThreeVector Precoil = ptrk;

  // Get IVRate class
  G4CMPInterValleyRate* ivpro = dynamic_cast<G4CMPInterValleyRate*>(GetRateModel());

  // Make sure it's not empty
  if (ivpro == nullptr) { 
    G4cout << "NO G4CMPInterValleyRate* " << G4endl;
    return &aParticleChange;
  }

  // Get IV rates 
  else if (ivpro != nullptr) {
    std::vector<G4double> probabilities = ivpro->GetIVProb();
    std::vector<G4double> cumulatives(probabilities.size());
    G4double totalIVRate = 0.;

    // Calculate total IV rate
    for (auto& element : probabilities) { totalIVRate += element; }

    // Don't do anything if IV rate is 0 
    if (totalIVRate == 0) { 
      G4cout << "NO IV Rate " << G4endl;
      return &aParticleChange;
    }

    // Do IV scattering if IV rate is not 0 
    else if (totalIVRate !=0) {

      // Normalize IV rates 
      for (auto& element : probabilities) { element /= totalIVRate; }

      // Calculate weight of each rates to total rate
      cumulatives[0] = probabilities[0];
      for (size_t i = 1; i < probabilities.size(); i++) {
        cumulatives[i] = cumulatives[i - 1] + probabilities[i];
      }

      // Choose random rate index from weighted step distribution
      G4double ivrandom = G4UniformRand();
      G4int Ephononi = 0;
      for (size_t i = 0; i < cumulatives.size(); i++) {
        if (ivrandom < cumulatives[i]) {
          Ephononi = i;
          break;
        }
      }

      // Get track's kinematics in current valley
      G4double Etrk = GetKineticEnergy(aTrack);
      G4double Ephonon = theLattice->GetIVEnergy(Ephononi);

      // Don't do anything if phonon energy > electron energy
      if (Ephonon > Etrk) {
        Ephonon=0;
        G4cout << "Ephonon = 0 " << G4endl;
        return &aParticleChange;
      }

      // Do IV scattering if phonon energy < electron energy
      if (Ephonon!=0) {

        // track's kinematics before IV scattering in current valley
        G4ThreeVector ktrk = theLattice->MapPtoK(valley, ptrk);
        G4ThreeVector kHV = 
          theLattice->EllipsoidalToSphericalTranformation(valley, ktrk);
        G4double kmag = kHV.mag();

        // Initialize some of the track's kinematics after IV scattering
        G4double phi_phonon=0, q=0, Erecoil=0, costheta=0;
        G4ThreeVector qvec, k_recoilHV, k_recoil; 

        Erecoil = Etrk - Ephonon; // electron energy after phonon emission

        // Get angle of phonon emission (0th or 1st order)
        G4double ivorder = theLattice->GetIVOrder(Ephononi); 
        if (ivorder==0) costheta = MakePhononThetaIV0Order(Etrk,Ephonon);
        if (ivorder==1) costheta = MakePhononThetaIV1Order(Etrk,Ephonon);
        phi_phonon = G4UniformRand() * twopi;

        // Select randomly one of the phonon momentum branch   
        if (G4UniformRand() <0.5) {
          q = kmag * costheta - kmag * 
              sqrt(costheta * costheta - Ephonon / Etrk);}
        else {
          q = kmag * costheta + kmag * 
              sqrt(costheta * costheta - Ephonon / Etrk);}

        // Compute the phonon momentum in HV frame
        G4ThreeVector kdir = kHV.unit();
        qvec = q * kdir;
        qvec.rotate(kdir.orthogonal(), acos(costheta));
        qvec.rotate(kdir, phi_phonon);

        // Compute the scattered electron in the HV frame
        k_recoilHV = kHV - qvec;
        // Making sure energy and momentum are conserved
        k_recoil  = theLattice->SphericalToEllipsoidalTranformation(valley, k_recoilHV);
        Precoil = theLattice->MapKtoP(valley, k_recoil);
        Precoil = theLattice->MapEkintoP(valley, Precoil.unit(),Erecoil);
        k_recoil = theLattice->MapPtoK(valley, Precoil);
        k_recoilHV = theLattice->EllipsoidalToSphericalTranformation(valley, k_recoil);           

        // IV f or g-type scattering
        G4String ivfgscat = theLattice->GetIVFGScattering(Ephononi); 

        // Find anti-valley
        G4int antivalley = 0;

        if (valley%2==0) {antivalley=valley+1;}
        if (valley%2!=0) {antivalley=valley-1;}

        // Picking anti-valley as new valley if g-type scattering
        if (ivfgscat=="g")  {valley=antivalley;}  

        // Picking at random any valley other than anti-valley if f-type scattering
        else if (ivfgscat=="f")  { 
          G4int fvalley = valley;
          while ( (fvalley==valley) || (fvalley==antivalley) ) {
            fvalley=ChangeValley(fvalley);
          }
          valley=fvalley;
        }

        // Picking a new valley at random if not f or g-type scattering
        else {
          valley = ChangeValley(valley);
        }

        // Calculate track's kinematics after IV scattering in new valley
        k_recoil = theLattice->SphericalToEllipsoidalTranformation(valley, k_recoilHV); 

        // Make sure electron is in valley and not anti-valley
        G4ThreeVector k0 = theLattice->RotateFromValley(valley, G4ThreeVector(1,0,0));
        if (k_recoil*k0/k_recoil.mag()/k0.mag() <0) {
          k_recoil=k_recoil*-1;
        }

        // Convert quasi-momentum to transport momentum
        Precoil = theLattice->MapKtoP(valley, k_recoil);

        // phonon wavector
        qvec = ktrk - k_recoil;
        q = qvec.mag();

        // Create real phonon to be propagated, with chosen polarization 
        G4int phononmode=G4PhononPolarization::Long; // LA default phonon mode
        G4String phononbranch=
          theLattice->GetIVPhononMode(Ephononi); // Get phonon mode 

        // If phonon mode is TA -> FT or ST 
        if (phononbranch=="TA") {
          phononmode=G4CMP::ChoosePhononPolarization(
            0,theLattice->GetSTDOS(),theLattice->GetFTDOS());
        }

        // If phonon mode is LO or TO keep -> LA 

        G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(aTrack)->SetValleyIndex(valley);

        G4double weight =
          G4CMP::ChoosePhononWeight(G4CMPConfigManager::GetLukeSampling());
        if (weight > 0.) {
          G4Track* phonon = G4CMP::CreatePhonon(aTrack,
                                                phononmode,
                                                qvec, Ephonon,
                                                aTrack.GetGlobalTime(),
                                                aTrack.GetPosition());
          // Secondary's weight has to be multiplicative with its parent's
          phonon->SetWeight(aTrack.GetWeight() * weight);
          if (verboseLevel > 1) {
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
        }
        else {
          aParticleChange.ProposeNonIonizingEnergyDeposit(Ephonon);
        }
      }
    }
  }

  // Rotate transport momentum to global frame
  RotateToGlobalDirection(Precoil);

  // Adjust track kinematics for new valley
  FillParticleChange(valley, Precoil);
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
