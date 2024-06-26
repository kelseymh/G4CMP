/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRecombinationRate.cc
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
//
#include "G4CMPBogoliubovQPRecombinationRate.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include <vector>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Recombination rate is computed using energy and the G4SCUtils class, upon which this is based
G4double G4CMPBogoliubovQPRecombinationRate::Rate(const G4Track& aTrack) const
{
  G4cout << "REL in BogoliubovQPRecombination rate Rate function." << G4endl;
  //Compute tau for recombination, and invert for rate
  G4double energy = GetKineticEnergy(aTrack);
  G4cout << "REL HereA in BogoliubovRecombination Rate" << G4endl;
  G4double tau_recombination = fTau0_qp*(this->GetTauAsAFunctionOfEnergy(fCurrentNormalizedTauRecombinationVsEnergy,"BogoliubovQP",energy));
  G4cout << "REL HereB in BogoliubovRecombination Rate" << G4endl;
  return (1.0/tau_recombination);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// If we arrive in a new lattice, either compute the recombination and then update the lookup table, or update
// the lookup table using a computed rate curve that already exists in the map
void G4CMPBogoliubovQPRecombinationRate::UpdateLookupTable(const G4LatticePhysical * theLat)
{
  //1. If the lattice doesn't exist in the lattice container associated with this process yet,
  //   add it and do the full calculation of the curves we care about, storing them in a map
  if( fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.count(theLat) == 0 ){
    G4cout << "Computing new lookup table for recombination process." << G4endl;
    fMap_physicalLattice_NormalizedTauRecombinationVsEnergy.emplace(theLat,ComputeNormalizedTauRecombinationVsEnergy());
    fCurrentNormalizedTauRecombinationVsEnergy = fMap_physicalLattice_NormalizedTauRecombinationVsEnergy[theLat];
  }    
  //2. If it does, just make the "active" functions the ones in the map
  else{  fCurrentNormalizedTauRecombinationVsEnergy = fMap_physicalLattice_NormalizedTauRecombinationVsEnergy[theLat]; }  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for recombination vs phonon energy
std::vector<std::vector<G4double> > G4CMPBogoliubovQPRecombinationRate::ComputeNormalizedTauRecombinationVsEnergy()
{
  std::vector<std::vector<G4double> > output;
  G4double deltaQPEnergyDivGap = (fMaxQPEnergyDivGap - fMinQPEnergyDivGap) / ((double)fQPEnergyBins);

  //Loop over all QP energy bins, and create a normalized tau for each of them
  for( int iB = 0; iB < fQPEnergyBins; ++iB ){
    
    //Define the qp energy for this bin    
    G4double qpEnergy = fGapEnergy * (deltaQPEnergyDivGap * (iB+0.5) + fMinQPEnergyDivGap);
    G4double prefactor = 2*pi / hbar_Planck / (1-FermiFactor(qpEnergy,fTeff));

    //Now do an integral over Omega, phonon energy
    int nW = 1000000;
    double minOmega = qpEnergy + fGapEnergy; //Minimum you can get is QP energy + another QP at gap edge
    double maxOmega = 1 * CLHEP::eV; // Go up to 1 eV
    double deltaOmega = (maxOmega - minOmega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++iW ){
      double Omega = minOmega + (iW + 0.5) * deltaOmega;
      double alpha2F = Omega*Omega;
      double energyTerm1 = (Omega-qpEnergy)/pow(pow(Omega-qpEnergy,2)-fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 = (1 + fGapEnergy*fGapEnergy / qpEnergy / (Omega-qpEnergy));
      double integrand = alpha2F * energyTerm1 * energyTerm2 * (BoseFactor(Omega,fTeff)+1) * (FermiFactor(Omega-qpEnergy,fTeff));
      integral += integrand * deltaOmega;
    }
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 2 / pi / pow(k_Boltzmann*fTcrit,3);
    double normalizedTau = 1.0 / inverseTau / tau0;

    //For now we use this. But can optimize by making this of an array instead of an std::vector
    std::vector<G4double> element;
    element.push_back(qpEnergy);
    element.push_back(normalizedTau);
    output.push_back(element);
  }

  //This is only for debugging, and is temporary.
  SaveBogoliubovRecombinationTauVsPhononEnergyToLogFile(output);

  
  return output;
}


//REL NEED TO SWAP DIRECT CALLS TO PROTECTED DATA MEMBERS WITH CONST ACCESS FUNCTIONS

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for pairbreaking vs phonon energy
void G4CMPBogoliubovQPRecombinationRate::SaveBogoliubovRecombinationTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc)
{
  std::ofstream outfile;
  outfile.open("/Users/ryanlinehan/QSC/Sims/Geant4/scRebuild-build/BogoliubovRecombinationTauVsEnergy.txt");
  for( int iE = 0; iE < theFunc.size(); ++iE ){
    outfile << theFunc[iE][0] << " " << theFunc[iE][1] << std::endl;
  }
  outfile.close();
  return;
}
