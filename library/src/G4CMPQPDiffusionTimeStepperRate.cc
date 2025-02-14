/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPBogoliubovQPRecombinationRate.cc
/// \brief Output a rate for an artificial time step for QPs undergoing diffusion
//
#include "G4CMPQPDiffusionTimeStepperRate.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include <vector>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4CMPQPDiffusionTimeStepperRate::Rate(const G4Track& aTrack) const
{
  G4cout << "REL in QPDiffusionTimeStepper rate Rate function." << G4endl;
  //Compute tau for recombination, and invert for rate
  double tau_nextScatter = 100000000 * CLHEP::ns; //Temporarily hardcoded
  return (1.0/tau_nextScatter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4CMPQPDiffusionTimeStepperRate::UpdateLookupTable(const G4LatticePhysical * theLat)
{
  G4cout << "In the updateLookupTable function, QPDiffusionTimeStepper. Nothing happens here because this class does not (yet) need a lookup table (but this function is generic to G4CMPVProcesses." << G4endl;
  return;
}
