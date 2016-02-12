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
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackInformation.hh"
#include "G4PhononReflection.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononPolarization.hh"
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
#include "Randomize.hh"


G4PhononReflection::G4PhononReflection(const G4String& aName)
  : G4VPhononProcess(aName, fPhononReflection),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4VParticleChange* G4PhononReflection::PostStepDoIt(const G4Track& aTrack,
                                                    const G4Step& aStep) {
  // This process handles the interaction of phonons with
  // boundaries. Implementation of this class is highly geometry
  // dependent.Currently, phonons are killed when they reach a
  // boundary. If the other side of the boundary was Al, a hit is
  // registered.
  aParticleChange.Initialize(aTrack);
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  G4StepPoint* preStepPoint = aStep.GetPreStepPoint();

  // Do nothing but return if the step is not boundary limited, or if
  // the step is too short. This is to allow actual reflection where after
  // the first boundary crossing a second, infinitesimal, step occurs crossing
  // back into the original volume
  if (postStepPoint->GetStepStatus() != fGeomBoundary ||
      aTrack.GetStepLength() <= kCarTolerance/2.) {
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
  // FIXME: Shouldn't this test be redundant with the kCarTolerance test?
  G4LatticePhysical*
    volLattice = G4LatticeManager::GetLatticeManager()->GetLattice(thePrePV);
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
  G4CMPSurfaceProperty* borderSurface = nullptr;
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

  if (!borderSurface) {
    if (verboseLevel>1) {
      G4cerr << GetProcessName() << ": No surface properties defined for "
             << thePrePV->GetName() << " to "  << thePostPV->GetName()
             << G4endl;
    }

    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  G4ThreeVector surfNorm =
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
  G4CMPTrackInformation* trackInfo = static_cast<G4CMPTrackInformation*>(
    aTrack.GetAuxiliaryTrackInformation(fPhysicsModelID));
  G4int nRefl = trackInfo->GetReflectionCount();
  if (maxRefl<0 || nRefl < maxRefl) {
    trackInfo->IncrementReflectionCount();
    return DoReflection(aStep, surfNorm, borderSurface, trackInfo);
  } else {
    if (verboseLevel) {
      G4cout << GetProcessName() << " WARNING: Phonon has reflected "
             << maxRefl << " times. Track being killed." << G4endl;
    }
    G4double eKin = aTrack.GetKineticEnergy();
    aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  return &aParticleChange; 
}

G4double G4PhononReflection::GetMeanFreePath(const G4Track&aTrack,
                                             G4double prevStepLength,
                                             G4ForceCondition* condition) {
  // Always return DBL_MAX and Forced. This ensures that the process is
  // called at the end of every step. In PostStepDoIt the process
  // decides whether the step encountered a volume boundary and a
  // reflection should be applied
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange*
G4PhononReflection::DoReflection(const G4Step& aStep,
                                 const G4ThreeVector& surfNorm,
                                 const G4CMPSurfaceProperty* surfProp,
                                 G4CMPTrackInformation* trackInfo) {
  if (verboseLevel>1)
    G4cout << GetProcessName() << ": Track reflected" << G4endl;

  G4ThreeVector k = trackInfo->GetPhononK();
  if (verboseLevel>2)
    G4cout << " Old momentum direction " << k.unit() << G4endl;

  // Grab phonon properties from the surface
  G4MaterialPropertiesTable*
    phonPropTable = const_cast<G4MaterialPropertiesTable*>(
      surfProp->GetPhononMaterialPropertiesTablePointer());
  G4double absProb       = phonPropTable->GetConstProperty("absProb");
  G4double specProb      = phonPropTable->GetConstProperty("specProb");

  if (G4UniformRand() < absProb) {
    // Kill the track, but potentially make new reflected phonons.
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    G4double Ekin = aStep.GetPostStepPoint()->GetKineticEnergy();
    std::vector<G4double> phononEnergies = KaplanPhononQP(Ekin, phonPropTable);

    for (G4double E : phononEnergies) {
      // TODO: Map E to K properly.
      G4double kmag = k.mag()*E/Ekin;
      G4ThreeVector reflectedKDir = k.unit();
      G4int pol = ChoosePolarization();
      do {
        reflectedKDir = LambertRotation(surfNorm);
      } while (!ReflectionIsGood(pol, kmag*reflectedKDir, surfNorm));

      CreatePhonon(pol,kmag*reflectedKDir,E);
    }
  } else {
    G4ThreeVector reflectedKDir = k.unit();
    G4int pol = GetPolarization(aStep.GetTrack());
    do {
      if (G4UniformRand() < specProb) {
        // Specular reflecton reverses momentum along normal
        G4double momNorm = reflectedKDir * surfNorm;
        reflectedKDir -= 2. * momNorm * surfNorm;
      } else {
        reflectedKDir = LambertRotation(surfNorm);
      }
      k = k.mag()*reflectedKDir;
    } while (!ReflectionIsGood(pol, k, surfNorm));

    if (verboseLevel>2)
      G4cout << " New momentum direction " << reflectedKDir << G4endl;

    G4ThreeVector vdir = theLattice->MapKtoVDir(pol, k);
    G4double v = theLattice->MapKtoV(pol, k);
    trackInfo->SetPhononK(k);
    aParticleChange.ProposeVelocity(v);
    aParticleChange.ProposeMomentumDirection(vdir);
  } // if absorbed

  return &aParticleChange;
}

G4ThreeVector
G4PhononReflection::LambertRotation(const G4ThreeVector& surfNorm) {
  G4double phi = 2.0*pi*G4UniformRand();
  G4double theta = acos(2.0*G4UniformRand() - 1.0) / 2.0;

  G4ThreeVector refl = -surfNorm;
  refl = refl.rotate(surfNorm.orthogonal(), theta);
  refl = refl.rotate(surfNorm, phi);
  return refl;
}

G4bool G4PhononReflection::ReflectionIsGood(G4int polarization,
                                            const G4ThreeVector& k,
                                            const G4ThreeVector& surfNorm) {
  G4ThreeVector vDir = theLattice->MapKtoVDir(polarization, k);
  return vDir.dot(surfNorm) < 0.0;
}

std::vector<G4double> G4PhononReflection::KaplanPhononQP(G4double energy,
                                        G4MaterialPropertiesTable* prop) {
  std::vector<G4double> reflectedEnergies;
  G4double gapEnergy     = prop->GetConstProperty("gapEnergy");
  G4double lowQPLimit    = prop->GetConstProperty("lowQPLimit");

  //For the phonon to not break a Cooper pair, it must go 2*thickness
  G4double frac = 2.0;

  G4double phononEscapeProb = CalcEscapeProbability(energy, frac, prop);

  G4double EDep = 0.;
  if (energy > 2.0*gapEnergy && G4UniformRand() > phononEscapeProb) {
    std::vector<G4double> qpEnergies;
    std::vector<G4double> phonEnergies{energy};
    while (qpEnergies.size()>0 || phonEnergies.size()>0) {
      if (phonEnergies.size()>0) {
        // Partition the phonons' energies into quasi-particles according to
        // a PDF defined in CalcQPEnergies().
        // NOTE: Both energy vectors mutate.
        EDep += CalcQPEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (qpEnergies.size()>0) {
        // Quasiparticles can also excite phonons.
        // NOTE: Both energy vectors mutate.
        EDep += CalcPhononEnergies(gapEnergy, lowQPLimit, phonEnergies, qpEnergies);
      }
      if (phonEnergies.size()>0) {
        // Some phonons will escape back into the crystal.
        // NOTE: Both energy vectors mutate.
        CalcReflectedPhononEnergies(prop, phonEnergies, reflectedEnergies);
      }
    }
  } else {
    reflectedEnergies.push_back(energy);
  }

  aParticleChange.ProposeNonIonizingEnergyDeposit(EDep);
  return reflectedEnergies;
}

G4double G4PhononReflection::CalcEscapeProbability(G4double energy,
                                        G4double thicknessFrac,
                                        G4MaterialPropertiesTable* prop) {
  G4double gapEnergy = prop->GetConstProperty("gapEnergy");
  G4double phononLifetime = prop->GetConstProperty("phononLifetime");
  G4double phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
  G4double vSound = prop->GetConstProperty("vSound");
  G4double thickness = prop->GetConstProperty("filmThickness");

  G4double mfp = vSound * phononLifetime / (1. + phononLifetimeSlope * (
    energy/gapEnergy - 2.));
  return exp(-2.* thicknessFrac * thickness/mfp);
}

G4double G4PhononReflection::CalcQPEnergies(G4double gapEnergy,
                                            G4double lowQPLimit,
                                            std::vector<G4double>& phonEnergies,
                                            std::vector<G4double>& qpEnergies) {
  // Each phonon gives all of its energy to the qp pair it breaks.
  G4double EDep = 0.;
  for (G4double E: phonEnergies) {
    G4double qpE = QPEnergyRand(gapEnergy, E);
    if (qpE >= lowQPLimit*gapEnergy)
      qpEnergies.push_back(qpE);
    else
      EDep += qpE;

    if (E-qpE >= lowQPLimit*gapEnergy)
      qpEnergies.push_back(E-qpE);
    else
      EDep += E-qpE;
  }

  phonEnergies.clear();
  return EDep;
}


G4double G4PhononReflection::CalcPhononEnergies(G4double gapEnergy,
                                            G4double lowQPLimit,
                                            std::vector<G4double>& phonEnergies,
                                            std::vector<G4double>& qpEnergies) {
  // NOTE: Phonons with low energy will not be seen by the detector, so we
  // don't record those energies and just "lose" those phonons.
  // Have a reference in for loop b/c qp doesn't give all of its energy away.
  G4double EDep = 0.;
  std::vector<G4double> newQPEnergies;
  for (G4double& E: qpEnergies) {
    // NOTE: E mutates in PhononEnergyRand.
    G4double phonE = PhononEnergyRand(gapEnergy, E);
    if (phonE >= 2.0*gapEnergy)
      phonEnergies.push_back(phonE);
    if (E >= lowQPLimit*gapEnergy)
      newQPEnergies.push_back(E);
    else
      EDep += E;
  }

  qpEnergies.swap(newQPEnergies);
  return EDep;
}

void G4PhononReflection::CalcReflectedPhononEnergies(
                                     G4MaterialPropertiesTable* prop,
                                     std::vector<G4double>& phonEnergies,
                                     std::vector<G4double>& reflectedEnergies) {
  // There is a 50% chance that a phonon is headed away from (toward) substrate
  std::vector<G4double> newPhonEnergies;
  for (G4double E : phonEnergies) {
    G4double frac = (G4UniformRand() < 0.5 ? 0.5 : 1.5);
    if (G4UniformRand() < CalcEscapeProbability(E, frac, prop))
      reflectedEnergies.push_back(E);
    else
      newPhonEnergies.push_back(E);
  }
  phonEnergies.swap(newPhonEnergies);
}

G4double G4PhononReflection::QPEnergyRand(G4double gapEnergy, G4double Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy - E') + gapEnergy*gapEnergy)
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy) *
  //                ((Energy - E')*(Energy - E') - gapEnergy*gapEnergy));
  // The shape of the PDF is like a U, so the max values are at the endpoints:
  // E' = gapEnergy and E' = Energy - gapEnergy

  // Add buffer so first/last bins don't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + (Energy-2.*gapEnergy)/BUFF;
  G4double xmax = gapEnergy + (Energy-2.*gapEnergy)*(BUFF-1.)/BUFF;

  G4double ymax = (xmin*(Energy-xmin) + gapEnergy*gapEnergy)
                  /
                  sqrt((xmin*xmin - gapEnergy*gapEnergy) *
                       ((Energy-xmin)*(Energy-xmin) - gapEnergy*gapEnergy));

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest) + gapEnergy*gapEnergy)
                  /
                  sqrt((xtest*xtest - gapEnergy*gapEnergy) *
                       ((Energy-xtest)*(Energy-xtest) - gapEnergy*gapEnergy))) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  return xtest;
/*
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // PDF is shaped like a capital U, so a rejection method would be slow.
  // Let's numerically calculate a CDF and throw a random to find E.

  const G4int BINS = 1000;
  G4double energyArr[BINS];
  G4double cdf[BINS];
  G4double pdfSum = 0.;
  for (size_t i = 0; i < BINS; ++i) {
    // Add 1 to i so first bin doesn't give zero denominator in pdfSum
    // Add 1 to BINS so last bin doesn't give zero denominator in pdfSum
    energyArr[i] = gapEnergy + (Energy-2.*gapEnergy) * (i+1.)/(BINS+1.);
    pdfSum += (energyArr[i]*(Energy - energyArr[i]) + gapEnergy*gapEnergy)
              /
              sqrt((energyArr[i]*energyArr[i] - gapEnergy*gapEnergy) *
                   ((Energy - energyArr[i])*(Energy - energyArr[i]) -
                    gapEnergy*gapEnergy));
    cdf[i] = pdfSum;
  }

  G4double u = G4UniformRand();

  size_t index = 0;
  for (; index < BINS; ++index) { //Combine normalization and search loops
    if (cdf[index]/pdfSum >= u)
      break;
  }

  return energyArr[index];
*/
}

G4double
G4PhononReflection::PhononEnergyRand(G4double gapEnergy, G4double& Energy) {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
  //           /
  //           sqrt((E'*E' - gapEnergy*gapEnergy);

  // Add buffer so first bin doesn't give zero denominator in pdfSum
  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + gapEnergy/BUFF;
  G4double xmax = Energy;

  G4double ymax = (xmin*(Energy-xmin)*(Energy-xmin) *
                    (xmin-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xmin*xmin - gapEnergy*gapEnergy);

  G4double ytest = G4UniformRand()*ymax;
  G4double xtest = G4UniformRand()*(xmax-xmin) + xmin;
  while (ytest > (xtest*(Energy-xtest)*(Energy-xtest) *
                    (xtest-gapEnergy*gapEnergy/Energy)) /
                  sqrt(xtest*xtest - gapEnergy*gapEnergy)) {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  }

  G4double phononE = Energy - xtest;
  Energy = xtest;
  return phononE;
}
