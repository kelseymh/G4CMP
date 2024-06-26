/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRecombinationRate.cc
/// \brief Compute rate for phonon breaking CPs into 2 QPs.
//
#include "G4CMPBogoliubovQPRadiatesPhononRate.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include <vector>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Recombination rate is computed using energy and the G4SCUtils class, upon which this is based
G4double G4CMPBogoliubovQPRadiatesPhononRate::Rate(const G4Track& aTrack) const
{
  G4cout << "REL in BogoliubovQPRadiatesPhonon rate Rate function." << G4endl;
  //Compute tau for recombination, and invert for rate
  G4double energy = GetKineticEnergy(aTrack);
  G4cout << "REL HereA in BogoliubovRadiatesPhonon Rate" << G4endl;
  G4double tau_scattering = fTau0_qp*(this->GetTauAsAFunctionOfEnergy(fCurrentNormalizedTauQPRadiatesPhononVsEnergy,"BogoliubovQP",energy));
  G4cout << "REL HereB in BogoliubovRecombination Rate" << G4endl;
  return (1.0/tau_scattering);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// If we arrive in a new lattice, either compute the scattering rate and then update the lookup table, or update
// the lookup table using a computed rate curve that already exists in the map
void G4CMPBogoliubovQPRadiatesPhononRate::UpdateLookupTable(const G4LatticePhysical * theLat)
{
  //1. If the lattice doesn't exist in the lattice container associated with this process yet,
  //   add it and do the full calculation of the curves we care about, storing them in a map
  if( fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy.count(theLat) == 0 ){
    G4cout << "Computing new lookup table for phonon radiation process." << G4endl;
    fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy.emplace(theLat,ComputeNormalizedTauQPRadiatesPhononVsEnergy());
    fCurrentNormalizedTauQPRadiatesPhononVsEnergy = fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy[theLat];
  }    
  //2. If it does, just make the "active" functions the ones in the map
  else{  fCurrentNormalizedTauQPRadiatesPhononVsEnergy = fMap_physicalLattice_NormalizedTauQPRadiatesPhononVsEnergy[theLat]; }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for scattering vs phonon energy. Note that this is JUST for QPs radiating phonons
// at the moment
std::vector<std::vector<G4double> > G4CMPBogoliubovQPRadiatesPhononRate::ComputeNormalizedTauQPRadiatesPhononVsEnergy()
{  
  std::vector<std::vector<G4double> > output;
  G4double deltaQPEnergyDivGap = (fMaxQPEnergyDivGap - fMinQPEnergyDivGap) / ((double)fQPEnergyBins);

  //Loop over all QP energy bins, and create a normalized tau for each of them
  for( int iB = 0; iB < fQPEnergyBins; ++iB ){

    //Define the qp energy for this bin    
    G4double qpEnergy = fGapEnergy * (deltaQPEnergyDivGap * (iB+0.5) + fMinQPEnergyDivGap);
    G4double prefactor = 2*pi / hbar_Planck / (1-FermiFactor(qpEnergy,fTeff)); //Removing Z0 since it's canceled later
    
    //Now do an integral over Omega, phonon energy
    int nW = 1000; //Default is 10k
    G4double minOmega = 0;
    G4double maxOmega = qpEnergy - fGapEnergy;
    G4double deltaOmega = (maxOmega - minOmega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++iW ){
      double Omega = minOmega + (iW + 0.5) * deltaOmega;
      double alpha2F = Omega*Omega; // Omitting b since it's divided out later
      double energyTerm1 = (qpEnergy-Omega)/pow(pow(qpEnergy-Omega,2)-fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 = (1 - fGapEnergy*fGapEnergy / qpEnergy / ( qpEnergy - Omega ));
      double integrand = alpha2F * energyTerm1 * energyTerm2 * (BoseFactor(Omega,fTeff) + 1) * (1 - FermiFactor(qpEnergy-Omega,fTeff));
      integral += integrand * deltaOmega;
    }
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 2 / pi / pow(k_Boltzmann*fTcrit,3); //Again omitting b and Z0 since this is where they're divided out
    double normalizedTau = 1.0/inverseTau / tau0;  
    
    //For now we use this. But can optimize by making this of an array instead of an std::vector
    std::vector<G4double> element;
    element.push_back(qpEnergy);
    element.push_back(normalizedTau);
    output.push_back(element);
  }

  //This is only for debugging, and is temporary.
  SaveBogoliubovQPRadiatesPhononTauVsPhononEnergyToLogFile(output);

  
  return output;
}


//REL NEED TO SWAP DIRECT CALLS TO PROTECTED DATA MEMBERS WITH CONST ACCESS FUNCTIONS

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for pairbreaking vs phonon energy
void G4CMPBogoliubovQPRadiatesPhononRate::SaveBogoliubovQPRadiatesPhononTauVsPhononEnergyToLogFile(std::vector<std::vector<G4double> > theFunc)
{
  std::ofstream outfile;
  outfile.open("/Users/ryanlinehan/QSC/Sims/Geant4/scRebuild-build/BogoliubovQPRadiatesPhononTauVsEnergy.txt");
  for( int iE = 0; iE < theFunc.size(); ++iE ){
    outfile << theFunc[iE][0] << " " << theFunc[iE][1] << std::endl;
  }
  outfile.close();
  return;
}
