/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPQPRecombinationRate.cc
/// \brief Compute rate for QP recombining with ambient bath QP into 2Delta
///  phonon
//
// 20260826  G4CMP-662: Splitting off functionality to check rate lookup table
//       for a given input lattice


#include "G4CMPQPRecombinationRate.hh"
#include "G4LatticePhysical.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4CMPSCUtils.hh"
#include <vector>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4CMPQPRecombinationRate::Rate(const G4Track& aTrack) const {
  return Rate(aTrack,0);
}

// Recombination rate is computed using energy and the G4SCUtils class, upon
// which this is based. If there is an additional lattice provided, then compute
// the rate based on that lattice
G4double G4CMPQPRecombinationRate::
Rate(const G4Track& aTrack, const G4LatticePhysical* theLat) const {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationRate::Rate --" << G4endl;
  }
  
  //First, check to see if the input lattice is zero. If it is, then we're
  //going to use the "current" lattice to do the rate calculation.
  bool useInputLat = false;
  if (theLat != 0) useInputLat = true;

  //Put checks to see if parameters are defined HERE -- this happens before
  //the calls to the vector but has access to tau0_qp, etc.
  if (useInputLat) {
    if (!CheckToSeeSCParametersSet(theLat)) return 0;
  } else {
    if (!CheckToSeeSCParametersSet()) return 0;
  }

  //Boolean for checking to see if we're trying to access below our minimum
  //energy (in the case of a turnaround step)
  bool thisEnergyBelowUsableRange = false;
  
  //Get the energy of the quasiparticle
  G4double energy = GetKineticEnergy(aTrack);

  //Compute the tau using either the current lattice's information or the
  //lattice information of the lattice provided as an argument.
  G4double tau_recombination = DBL_MAX;
  if (!useInputLat) {
    tau_recombination = fTau0_qp*
      (this->GetTauAsAFunctionOfEnergy
       (fCurrentNormalizedTauRecombinationVsEnergy,
        "QP",energy,thisEnergyBelowUsableRange));
  }
  else {

    //As a sanity check, make SURE that we have this lattice in our map of
    //vectors that give tau vs energy
    if (fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.count(theLat)
        == 0) {
      G4ExceptionDescription msg;
      msg << "For the input lattice, " << theLat->GetLattice()->GetName()
          << ", we noticed that we don't have an already-built normalized tau"
          << " versus energy. This is problematic. Fix.";
      G4Exception("G4CMPQPRecombinationRate::Rate",
                  "QPRecombinationRate000",FatalException, msg);
      return 0;
    }
    tau_recombination = theLat->GetSCTau0qp()*
      (this->GetTauAsAFunctionOfEnergy
       (fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.at(theLat),
        "QP",energy,thisEnergyBelowUsableRange,theLat));
  }

  //If we're below the usable range
  if (thisEnergyBelowUsableRange) {

    //Debugging
    if (verboseLevel > 5) {
      G4cout << "R Function Point A | In Rate calculation for QPRecombination,"
             << "this energy " << energy << " is below the usable range. "
             << "Returning a zero rate." << G4endl;
    }
    return 0;
  }

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "R Function Point B | tau_recombination: " << tau_recombination
           << G4endl;
    G4cout << "R Function Point B | returning 1.0/tau_recombination."
           << G4endl;
  }
  return (1.0/tau_recombination);
}

//This is meant to ensure that when we attempt to calculate a rate, we actually
//have the correct parameters set for this material, so that we exercise some
//control over the rate calculation.
bool G4CMPQPRecombinationRate::
CheckToSeeSCParametersSet(const G4LatticePhysical * theLat) const {

  //If a lattice is provided, we want to ask about that lattice info. If one is
  //not, we want to ask about the lattice info currently loaded in
  if (theLat == 0) {

    //Check for the gap0energy, Tcrit, Teff, and Tau0qp. If all of these aren't
    //set, return false. However, if any subset of them are set, then throw a
    //flag--means that someone may just have forgot one of them.
    if (fGap0Energy==0 || fTau0_qp == DBL_MAX || fTcrit == 0 || fTeff == 0) {

      //Means the whole material likely not set -- this is sometimes expected
      //during normal operation, so don't worry too much here.
      if (fGap0Energy==0 && fTau0_qp == DBL_MAX && fTcrit == 0 && fTeff == 0) {
        return false;
      } else {
        //^Means that the material is partially set -- this is probably a
        //mistake
        G4ExceptionDescription msg;
        msg << "Noticed that in the rate calculation step for the QP "
            << "recombination process, you have incorrectly defined or omitted "
            << "the Gap0Energy parameter, the Tcrit parameter, the Teff "
            << "parameter, or the Tau0qp parameter. In other words, you don't "
            << "have enough input information in your config.txt file to run "
            << "the recombination physics correctly.";
        G4Exception("G4CMPQPRecombinationRate::CheckToSeeSCParametersSet",
                    "QPRecombinationRate001",JustWarning, msg);
        return false;
      }
    }
    return true;
  }
  else {
    //Check for the gap0energy, Tcrit, Teff, and Tau0qp. If all of these aren't
    //set, return false. However, if any subset of them are set, then throw a
    //flag--means that someone may just have forgot one of them.
    if (theLat->GetSCDelta0()==0 || theLat->GetSCTau0qp() == DBL_MAX ||
        theLat->GetSCTcrit() == 0 || theLat->GetSCTeff() == 0) {

      //Means the whole material likely not set -- this is sometimes expected
      //during normal operation, so don't worry too much here.
      if (theLat->GetSCDelta0()==0 && theLat->GetSCTau0qp() == DBL_MAX &&
          theLat->GetSCTcrit() == 0 && theLat->GetSCTeff() == 0) {
        return false;
      } else {
        //^Means that the material is partially set -- this is probably a
        //mistake
        G4ExceptionDescription msg;
        msg << "Noticed that in the rate calculation step for the QP "
            << "recombination process, you have incorrectly defined or omitted "
            << "the Gap0Energy parameter, the Tcrit parameter, the Teff "
            << "parameter, or the Tau0qp parameter. In other words, you don't "
            << "have enough input information in your config.txt file to run "
            << "the recombination physics correctly.";
        G4Exception("G4CMPQPRecombinationRate::CheckToSeeSCParametersSet",
                    "QPRecombinationRate002",JustWarning, msg);
        return false;
      }
    }
    return true;
  }
}

// If we arrive in a new lattice, either compute the recombination and then
// update the lookup table, or update the lookup table using a computed rate
// curve that already exists in the map
void
G4CMPQPRecombinationRate::UpdateLookupTable(const G4LatticePhysical * theLat) {

  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationRate::UpdateLookupTable --" << G4endl;
    G4cout << "ULT Function Point A | QP recombination." << G4endl;
  } 
  
  //1. If the lattice doesn't exist in the lattice container associated with
  //   this process yet, add it and do the full calculation of the curves we
  //   care about, storing them in a map
  if (!CheckLookupTableForLat(theLat)) {
    if (verboseLevel > 5) {
      G4cout << "Computing new lookup table for recombination process, lattice "
             << "name: " << theLat->GetLattice()->GetName() << G4endl;
    }
    fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.
      emplace(theLat,ComputeNormalizedTauRecombinationVsEnergy());
    fCurrentNormalizedTauRecombinationVsEnergy =
      fMap_physicalLattice_NormalizedTauRecombinationVsEnergy[theLat];
  } else {
    //2. ^If it does, just make the "active" functions the ones in the map    
    fCurrentNormalizedTauRecombinationVsEnergy =
      fMap_physicalLattice_NormalizedTauRecombinationVsEnergy[theLat];
  }  
}

bool G4CMPQPRecombinationRate::
CheckLookupTableForLat(const G4LatticePhysical * theLat) const {
  //1. If the lattice doesn't exist in the lattice container associated with
  //   this process yet, throw a flag
  if(fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.count(theLat)
     == 0){ return false; }
  return true;
}


// Construct the lookup table for normalized tau for recombination vs phonon
// energy
std::vector<std::vector<G4double> >
G4CMPQPRecombinationRate::ComputeNormalizedTauRecombinationVsEnergy() {
  
  //Debugging
  if (verboseLevel > 5) {
    G4cout << "-- G4CMPQPRecombinationRate::"
           << "ComputeNormalizedTauRecombinationVsEnergy --" << G4endl;
    G4cout << "CNTRVE Function Point A | In the calculation of a normalized "
           << "tauQP vs energyQP, recombination." << G4endl;
  } 
  
  std::vector<std::vector<G4double> > output;
  G4double deltaQPEnergyDivGap =
    (fMaxQPEnergyDivGap - fMinQPEnergyDivGap) / ((double)fQPEnergyBins);
  
  //Loop over all QP energy bins, and create a normalized tau for each of them
  for (int iB = 0; iB < fQPEnergyBins; ++iB) {
    
    //Define the qp energy for this bin    
    G4double qpEnergy =
      fGapEnergy * (deltaQPEnergyDivGap * (iB+0.5) + fMinQPEnergyDivGap);
    G4double prefactor =
      2*pi / hbar_Planck / (1-FermiFactor(qpEnergy,fTeff));

    //Now do an integral over Omega, phonon energy
    //earlier was 1000000, but this was taking a loooong time
    int nW = 10000; 

    //Minimum you can get is QP energy + another QP at gap edge
    double minOmega = qpEnergy + fGapEnergy;
    
    //Somewhat arbitrary, but just "high" -- is supposed to represent infinity
    double maxOmega = fGapEnergy * fMaxQPEnergyDivGap * 10;    
    double deltaOmega = (maxOmega - minOmega) / ((double)nW);
    double integral = 0;
    for (int iW = 0; iW < nW; ++iW) {
      double Omega = minOmega + (iW + 0.5) * deltaOmega;
      double alpha2F = Omega*Omega;
      double energyTerm1 =
        (Omega-qpEnergy)/pow(pow(Omega-qpEnergy,2)-fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 =
        (1 + fGapEnergy*fGapEnergy / qpEnergy / (Omega-qpEnergy));
      double integrand =
        alpha2F * energyTerm1 * energyTerm2 * (BoseFactor(Omega,fTeff)+1) *
        (FermiFactor(Omega-qpEnergy,fTeff));
      integral += integrand * deltaOmega;
    }
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 2 / pi / pow(k_Boltzmann*fTcrit,3);
    double normalizedTau = 1.0 / inverseTau / tau0;

    //For now we use this. But can optimize by making this of an array instead
    //of an std::vector
    std::vector<G4double> element;
    element.push_back(qpEnergy);
    element.push_back(normalizedTau);
    output.push_back(element);
  }
  return output;
}


 
