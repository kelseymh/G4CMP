/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPSCPairBreakingRate.cc
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
//
#include "G4CMPSCPairBreakingRate.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4LatticePhysical.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Pairbreaking rate is computed using energy and the G4SCUtils class, upon which this is based
G4double G4CMPSCPairBreakingRate::Rate(const G4Track& aTrack) const
{
  G4cout << "REL HereA_SCPairBreakingRate" << G4endl;
  G4cout << "REL SCDelta for this lattice: " << (this->GetLattice())->GetSCDelta0() << G4endl;
  if( !CheckToSeeSCParametersSet() ){ return 0; }
  
  //Boolean for checking to see if we're trying to access below our minimum energy (in the case of a turnaround step)
  bool thisEnergyBelowUsableRange = false;
  
  //Need to check to see if the current lattice has a gap (i.e. IS a superconductor). If not, then we need return infinite for
  //this process, since it can in principle/code it can still run in the non-SC crystals.
  if( (this->GetLattice())->GetSCDelta0() <= 0.0 ){ return 0.0; }
  else{

    //Compute tau for pairbreaking, and invert for rate
    G4double energy = GetKineticEnergy(aTrack);
    G4cout << "REL HereAB_SCPairBreakingRate" << G4endl;
    G4double tau_pairbreaking = fTau0_ph*GetTauAsAFunctionOfEnergy(fCurrentNormalizedTauPairBreakingVsEnergy,"Phonon",energy,thisEnergyBelowUsableRange);
    G4cout << "REL HereB_SCPairBreakingRate" << G4endl;
    return (1.0/tau_pairbreaking);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//This is meant to ensure that when we attempt to calculate a rate, we actually have
//the correct parameters set for this material, so that we exercise some control over
//the rate calculation.
bool G4CMPSCPairBreakingRate::CheckToSeeSCParametersSet() const
{
  //Check for the gap0energy, Tcrit, Teff, and Tau0qp. If all of these aren't set, return false.
  //However, if any subset of them are set, then throw a flag--means that someone may just have forgot
  //one of them.
  if( fGap0Energy==0 || fTau0_ph == DBL_MAX || fTcrit == 0 || fTeff == 0 ){
    //Means the whole material likely not set -- this is sometimes expected during normal operation, so don't worry too much here.
    if( fGap0Energy==0 && fTau0_ph == DBL_MAX && fTcrit == 0 && fTeff == 0 ){
      return false;
    }
    //Means that the material is partially set -- this is probably a mistake
    else{
      G4ExceptionDescription msg;
      msg << "Noticed that in the rate calculation step for the SC Pairbreaking process, you have incorrectly defined or omitted the Gap0Energy parameter, the Tcrit parameter, the Teff parameter, or the Tau0ph parameter. In other words, you don't have enough input information in your config.txt file to run the pairbreaking physics correctly.";
      G4Exception("G4CMPSCPairbreakingRate::CheckToSeeSCParametersSet", "SCPairbreakingRate001",JustWarning, msg);
      return false;
    }
  }
  return true;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// If we arrive in a new lattice, either compute the pair breaking rate and then update the lookup table, or update
// the lookup table using a computed rate curve that already exists in the map
void G4CMPSCPairBreakingRate::UpdateLookupTable(const G4LatticePhysical * theLat)
{
  G4cout << "REL HereC_SCPairbreakingRate: updating lookup table." << G4endl;
  //1. If the lattice doesn't exist in the lattice container associated with this process yet,
  //   add it and do the full calculation of the curves we care about, storing them in a map
  if( fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy.count(theLat) == 0 ){
    G4cout << "Computing new lookup table for SC pairbreaking process." << G4endl;
    fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy.emplace(theLat,ComputeNormalizedTauPairBreakingVsEnergy());
    fCurrentNormalizedTauPairBreakingVsEnergy = fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy[theLat];
  }    
  //2. If it does, just make the "active" functions the ones in the map
  else{  fCurrentNormalizedTauPairBreakingVsEnergy = fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy[theLat]; }  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for pairbreaking vs phonon energy
std::vector<std::vector<G4double> > G4CMPSCPairBreakingRate::ComputeNormalizedTauPairBreakingVsEnergy()
{
  G4cout << "REL HereDA" << G4endl;
  std::vector<std::vector<G4double> > output;
  G4double deltaPhononEnergyDivGap = (fMaxPhononEnergyDivGap - fMinPhononEnergyDivGap) / ((double)fPhononEnergyBins);  
  
  //Loop over all phonon energy bins, and for each of them, compute the integral
  for( int iB = 0; iB < fPhononEnergyBins; ++iB ){

    //Define the phonon energy for this bin    
    G4double phononEnergy = fGapEnergy * (deltaPhononEnergyDivGap * (iB+0.5) + fMinPhononEnergyDivGap);

    //--------------------------------------------------
    //Now do the calculation for this phonon energy       
    //Initially, going to try this without specifying units -- shouuuuuld be fine but if there are problems this is the first place to debug.
    //Note that in what follows, omega indicates the QP energy integrated over.
    double prefactor = 4 * pi / hbar_Planck ; //Missing a factor of Alpha2 because it cancels later
    int nW = 20000;
    double minomega = fGapEnergy;
    double maxomega = phononEnergy - fGapEnergy;
    double deltaomega = (maxomega - minomega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++iW ){
      double omega = minomega + (iW+0.5)*deltaomega; //So that omega lands in the middle of a bin
      double energyTerm1 = 1.0 / pow(omega*omega - fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 = (omega*(phononEnergy-omega) + fGapEnergy*fGapEnergy) / pow( pow(phononEnergy-omega,2) - fGapEnergy*fGapEnergy, 0.5 );
      double integrand = energyTerm1 * energyTerm2 * (1-FermiFactor(phononEnergy-omega,fTeff) - FermiFactor(omega,fTeff));
      integral += integrand * deltaomega;
    }
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 4 / pi / pi / fGap0Energy; //Missing a factor of 1/Alpha2 because it cancels later
    double normalizedTau = 1.0 / inverseTau / tau0;

    //For now we use this. But can optimize by making this of an array instead of an std::vector
    std::vector<G4double> element;
    element.push_back(phononEnergy);
    element.push_back(normalizedTau);
    output.push_back(element);
  }
  G4cout << "REL HereDZ" << G4endl;


  //This is only for debugging, and is temporary.
  //  SavePairBreakingRateVsPhononEnergyToLogFile(output);

  
  return output;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for pairbreaking vs phonon energy
void G4CMPSCPairBreakingRate::SavePairBreakingRateVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc)
{
  std::ofstream outfile;
  outfile.open("/Users/ryanlinehan/QSC/Sims/Geant4/scRebuild-build/SCPairBreakingTauVsEnergy.txt");
  for( int iE = 0; iE < theFunc.size(); ++iE ){
    outfile << theFunc[iE][0] << " " << theFunc[iE][1] << std::endl;
  }
  outfile.close();
  return;
}



//REL NEED TO SWAP DIRECT CALLS TO PROTECTED DATA MEMBERS WITH CONST ACCESS FUNCTIONS
