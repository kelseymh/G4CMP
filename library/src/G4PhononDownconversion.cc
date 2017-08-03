/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4PhononDownconversion.cc
/// \brief Implementation of the G4PhononDownconversion class
//
// $Id$
//
// 20131111  Add verbose output for MFP calculation
// 20131115  Initialize data buffers in ctor
// 20140312  Follow name change CreateSecondary -> CreatePhonon
// 20140331  Add required process subtype code
// 20160624  Use GetTrackInfo() accessor
// 20161114  Use new PhononTrackInfo
// 20170620  Follow interface changes in G4CMPSecondaryUtils
// 20170802  Use G4CMP_DOWN_SAMPLE biasing with ChooseWeight(), move outside
//		of sub-functions.

#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4PhononDownconversion.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include <cmath>


G4PhononDownconversion::G4PhononDownconversion(const G4String& aName)
  : G4VPhononProcess(aName, fPhononDownconversion),
    fBeta(0.), fGamma(0.), fLambda(0.), fMu(0.) {
#ifdef G4CMP_DEBUG
  output.open("phonon_downsampling_stats", std::ios_base::app);
  if (output.good()) {
    output << "First Daughter Theta,Second Daughter Theta,First Daughter Energy [eV],Second Daughter Energy [eV],"
              "Decay Branch,First Daughter Weight,Second Daughter Weight,Parent Weight,"
              "Number of Outgoing Tracks,Parent Energy [eV]\n";
  } else {
    G4cerr << "Could not open phonon debugging output file!" << G4endl;}
#endif
  }

G4PhononDownconversion::~G4PhononDownconversion() {
#ifdef G4CMP_DEBUG
  output.close();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PhononDownconversion::GetMeanFreePath(const G4Track& aTrack,
						 G4double /*previousStepSize*/,
						 G4ForceCondition* condition) {
  //Determines mean free path for longitudinal phonons to split
  G4double A = theLattice->GetAnhDecConstant();
  G4double Eoverh = GetKineticEnergy(aTrack)/h_Planck;
  
  //Calculate mean free path for anh. decay
  G4double mfp = aTrack.GetVelocity()/(Eoverh*Eoverh*Eoverh*Eoverh*Eoverh*A);

  if (verboseLevel > 1)
    G4cout << "G4PhononDownconversion::GetMeanFreePath = " << mfp << G4endl;
  
  *condition = NotForced;
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange* G4PhononDownconversion::PostStepDoIt( const G4Track& aTrack,
							 const G4Step&) {
  aParticleChange.Initialize(aTrack);

  // Obtain dynamical constants from this volume's lattice
  fBeta   = theLattice->GetBeta() / (1e11*pascal);	// Make dimensionless
  fGamma  = theLattice->GetGamma() / (1e11*pascal);
  fLambda = theLattice->GetLambda() / (1e11*pascal);
  fMu     = theLattice->GetMu() / (1e11*pascal);

  //Destroy the parent phonon and create the daughter phonons.
  //74% chance that daughter phonons are both transverse
  //26% Transverse and Longitudinal
  if (G4UniformRand()>0.740) MakeLTSecondaries(aTrack);
  else MakeTTSecondaries(aTrack);

#ifdef G4CMP_DEBUG
  output << aTrack.GetWeight() << ','
         << aParticleChange.GetNumberOfSecondaries() << ','
         << aTrack.GetKineticEnergy()/eV << G4endl;
#endif

  // Only kill the track if downconversion actually happened
  if (aParticleChange.GetNumberOfSecondaries() > 0) {
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);    
  }

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PhononDownconversion::IsApplicable(const G4ParticleDefinition& aPD) {
  //Only L-phonons decay
  return (&aPD==G4PhononLong::PhononDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//probability density of energy distribution of L'-phonon in L->L'+T process

inline double G4PhononDownconversion::GetLTDecayProb(double d, double x) const {
  //d=delta= ratio of group velocities vl/vt and x is the fraction of energy in the longitudinal mode, i.e. x=EL'/EL
  return (1/(x*x))*(1-x*x)*(1-x*x)*((1+x)*(1+x)-d*d*((1-x)*(1-x)))*(1+x*x-d*d*(1-x)*(1-x))*(1+x*x-d*d*(1-x)*(1-x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//probability density of energy distribution of T-phonon in L->T+T process

inline double G4PhononDownconversion::GetTTDecayProb(double d, double x) const {  
  //dynamic constants from Tamura, PRL31, 1985
  G4double A = 0.5*(1-d*d)*(fBeta+fLambda+(1+d*d)*(fGamma+fMu));
  G4double B = fBeta+fLambda+2*d*d*(fGamma+fMu);
  G4double C = fBeta + fLambda + 2*(fGamma+fMu);
  G4double D = (1-d*d)*(2*fBeta+4*fGamma+fLambda+3*fMu);

  return (A+B*d*x-B*x*x)*(A+B*d*x-B*x*x)+(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)))*(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double G4PhononDownconversion::MakeLDeviation(double d, double x) const {
  //change in L'-phonon propagation direction after decay

  return std::acos((1+(x*x)-((d*d)*(1-x)*(1-x)))/(2*x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double G4PhononDownconversion::MakeTDeviation(double d, double x) const {
  //change in T-phonon propagation direction after decay (L->L+T process)
  
  return std::acos((1-x*x+d*d*(1-x)*(1-x))/(2*d*(1-x)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double G4PhononDownconversion::MakeTTDeviation(double d, double x) const {
  //change in T-phonon propagation direction after decay (L->T+T process)

  return std::acos((1-d*d*(1-x)*(1-x)+d*d*x*x)/(2*d*x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//Generate daughter phonons from L->T+T process
   
void G4PhononDownconversion::MakeTTSecondaries(const G4Track& aTrack) {
  //d is the velocity ratio vL/vT
  G4double d=1.6338;
  G4double upperBound=(1+(1/d))/2;
  G4double lowerBound=(1-(1/d))/2;

  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in first T phonon
  G4double x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
  G4double p = 1.5*G4UniformRand();
  while(p >= GetTTDecayProb(d, x*d)) {
    x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
    p = 1.5*G4UniformRand(); 
  }
  
  //using energy fraction x to calculate daughter phonon directions
  G4double theta1=MakeTTDeviation(d, x);
  G4double theta2=MakeTTDeviation(d, 1-x);
  G4ThreeVector dir1=G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack)->k();
  G4ThreeVector dir2=dir1;

  // FIXME:  These extra randoms change timing and causting outputs of example!
  //G4ThreeVector ran = G4RandomDirection();	// FIXME: Drop this line
  // Is this issue fixed by dropping the above line?
  
  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),theta1).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-theta2).rotate(dir2,ph);

  G4double E=GetKineticEnergy(aTrack);
  G4double Esec1 = x*E;
  G4double Esec2 = E-Esec1;

  // Make FT or ST phonons (0. means no longitudinal)
  G4int mode1 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
						theLattice->GetFTDOS());

  // Make FT or ST phonon (0. means no longitudinal)
  G4int mode2 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
						theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  // Always produce the secondaries.
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack.GetTouchable(), mode1,
				      dir1, Esec1, aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());
  G4Track* sec2 = G4CMP::CreatePhonon(aTrack.GetTouchable(), mode2,
                                      dir2, Esec2, aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());

  // Pick which secondary gets the weight randomly
#ifdef G4CMP_DEBUG
  output << theta1 << ',' << theta2 << ','
	 << sec1->GetKineticEnergy()/eV << ','
	 << sec2->GetKineticEnergy()/eV << ',';
#endif

  G4double bias = G4CMPConfigManager::GetDownconversionSampling();
  G4double weight = G4CMP::ChoosePhononWeight(bias);
  if (weight > 0.) {				// Produce both daughters
    aParticleChange.SetSecondaryWeightByProcess(true);
    sec1->SetWeight(aTrack.GetWeight()*weight); // Default weight
    sec2->SetWeight(aTrack.GetWeight()*weight);
#ifdef G4CMP_DEBUG
    output << "TT" << ',' << sec1->GetWeight() << ','
	   << sec2->GetWeight() << ',';
#endif

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec2);
    aParticleChange.AddSecondary(sec1);
  } else {			// Produce no daughters
    delete sec1;		// It would be better not to create/delete
    delete sec2;
#ifdef G4CMP_DEBUG
    output << "TT" << ',' << 0 << ',' << 0 << ',';
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//Generate daughter phonons from L->L'+T process
   
void G4PhononDownconversion::MakeLTSecondaries(const G4Track& aTrack) {
  //d is the velocity ratio vL/v
  G4double d=1.6338;
  G4double upperBound=1;
  G4double lowerBound=(d-1)/(d+1);
  
  /*
  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in L phonon
  G4double x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
  G4double p = 4.0*G4UniformRand();
  while(p >= GetLTDecayProb(d, x)) {
    x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
    p = 4.0*G4UniformRand(); 		     //4.0 is about the max in the PDF
  }
  */

  G4double u = G4UniformRand();
  G4double x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
  G4double q = 0;
  if (x <= upperBound && x >= lowerBound) q = 1/(upperBound-lowerBound);
  while (u >= GetLTDecayProb(d, x)/(2.8/(upperBound-lowerBound))) {
    u = G4UniformRand();
    x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
    if (x <= upperBound && x >= lowerBound) q = 1/(upperBound-lowerBound);
    else q = 0;
  }

  //using energy fraction x to calculate daughter phonon directions
  G4double thetaL=MakeLDeviation(d, x);
  G4double thetaT=MakeTDeviation(d, x);
  G4ThreeVector dir1=G4CMP::GetTrackInfo<G4CMPPhononTrackInfo>(aTrack)->k();
  G4ThreeVector dir2=dir1;

  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),thetaL).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-thetaT).rotate(dir2,ph);

  G4double E=GetKineticEnergy(aTrack);
  G4double Esec1 = x*E;
  G4double Esec2 = E-Esec1;

  // First secondary is longitudnal
  int mode1 = G4PhononPolarization::Long;

  // Make FT or ST phonon (0. means no longitudinal)
  G4int mode2 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
						theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  // Always produce the L mode phonon. Produce T mode phonon based on
  // biasing.
  G4Track* sec1 = G4CMP::CreatePhonon(aTrack.GetTouchable(), mode1,
				      dir1, Esec1, aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());
  G4Track* sec2 = G4CMP::CreatePhonon(aTrack.GetTouchable(), mode2,
                                      dir2, Esec2, aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());

#ifdef G4CMP_DEBUG
  output << thetaL << ',' << thetaT << ',' << sec1->GetKineticEnergy()/eV
	 << ',' << sec2->GetKineticEnergy()/eV << ',';
#endif

  G4double bias = G4CMPConfigManager::GetDownconversionSampling();
  G4double weight = G4CMP::ChoosePhononWeight(bias);
  if (weight > 0.) {				// Produce both daughters
    aParticleChange.SetSecondaryWeightByProcess(true);
    sec1->SetWeight(aTrack.GetWeight()/bias);
    sec2->SetWeight(aTrack.GetWeight()/bias);
#ifdef G4CMP_DEBUG
    output << "LT" << ',' << sec1->GetWeight() << ','
	   << sec2->GetWeight() << ',';
#endif

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec2);
    aParticleChange.AddSecondary(sec1);
  } else {			// Produce no daughters
    delete sec1;		// It would be better not to create/delete
    delete sec2;
#ifdef G4CMP_DEBUG
    output << "LT" << ',' << 0 << ',' << 0 << ',';
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

