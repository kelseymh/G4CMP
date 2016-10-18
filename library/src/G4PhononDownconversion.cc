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

#include "G4CMPTrackInformation.hh"
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
    fBeta(0.), fGamma(0.), fLambda(0.), fMu(0.) {;}

G4PhononDownconversion::~G4PhononDownconversion() {;}

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

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);    
       
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
  G4ThreeVector dir1=GetTrackInfo(aTrack)->GetPhononK();
  G4ThreeVector dir2=dir1;

  // FIXME:  These extra randoms change timing and causting outputs of example!
  //G4ThreeVector ran = G4RandomDirection();	// FIXME: Drop this line
  // Is this issue fixed by dropping the above line?
  
  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),theta1).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-theta2).rotate(dir2,ph);

  G4double E=GetKineticEnergy(aTrack);
  G4double Esec1 = x*E, Esec2 = E-Esec1;

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization1 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization2 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  // Always produce one of the secondaries. The other will be produced
  // based on track biasing values.
  G4Track* sec1 = CreatePhonon(polarization1, dir1, Esec1);
  G4Track* sec2 = CreatePhonon(polarization2, dir2, Esec2);

  // Pick which secondary gets the weight randomly
  if (G4UniformRand() < 0.5) {
    std::swap(sec1, sec2);
  }

  G4double weight1 = aTrack.GetWeight();
  G4double weight2 = weight1 * G4CMP::ChoosePhononWeight();
  if (weight2 > 0.) { // Produce both daughters
    aParticleChange.SetSecondaryWeightByProcess(true);
    sec1->SetWeight(weight1); // Default weight
    sec2->SetWeight(weight2);

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec2);
    aParticleChange.AddSecondary(sec1);
  } else { // Only produce one daughter
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(sec1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//Generate daughter phonons from L->L'+T process
   
void G4PhononDownconversion::MakeLTSecondaries(const G4Track& aTrack) {
  //d is the velocity ratio vL/v
  G4double d=1.6338;
  G4double upperBound=1;
  G4double lowerBound=(d-1)/(d+1);
  
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

  //using energy fraction x to calculate daughter phonon directions
  G4double thetaL=MakeLDeviation(d, x);
  G4double thetaT=MakeTDeviation(d, x);
  G4ThreeVector dir1=GetTrackInfo(aTrack)->GetPhononK();
  G4ThreeVector dir2=dir1;

  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),thetaL).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-thetaT).rotate(dir2,ph);

  G4double E=GetKineticEnergy(aTrack);
  G4double Esec1 = x*E, Esec2 = E-Esec1;

  // First secondary is longitudnal
  int polarization1 = G4PhononPolarization::Long;

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization2 = G4CMP::ChoosePhononPolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  // Always produce the L mode phonon. Produce T mode phonon based on
  // biasing.
  G4Track* sec1 = CreatePhonon(polarization1, dir1, Esec1);
  G4double weight1 = aTrack.GetWeight();
  G4double weight2 = weight1 * G4CMP::ChoosePhononWeight();
  if (weight2 > 0.) { // Produce both daughters
    G4Track* sec2 = CreatePhonon(polarization2, dir2, Esec2);

    aParticleChange.SetSecondaryWeightByProcess(true);
    sec1->SetWeight(weight1); // Default weight
    sec2->SetWeight(weight2);

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec2);
    aParticleChange.AddSecondary(sec1);
  } else { // Only produce L mode daughter
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(sec1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

