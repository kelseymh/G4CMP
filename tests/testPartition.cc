/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Usage: testPartition <Ehit> <Esample> <Lattice>
//
// Specify total hit energy (eV), downsampling threshold (eV),
// and lattice directory.  Geant4 material will be set as "G4_<Lattice>".
//
// NOTE: 10 keV energy deposit should produce ~3400 e/h pairs

#include "globals.hh"
#include "G4CMPEnergyPartition.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4Delete.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include <algorithm>
#include <stdlib.h>
#include <vector>

// Global variables for use in tests

namespace {
  G4CMPConfigManager* g4cmp = G4CMPConfigManager::Instance();
  G4LatticePhysical* lattice = 0;
  G4CMPEnergyPartition* partition = 0;
}

struct EPartInfo {
  G4double Ehit, Esamp, Echg, Ephon;
  G4int Nchg, Nphon;
};


// Test partitioning system one time

EPartInfo testPartition(G4double Ehit, G4double Esamp) {
  g4cmp->SetSamplingEnergy(Esamp);
  partition->DoPartition(Ehit, 0.);

  std::vector<G4PrimaryParticle*> prim;
  partition->GetPrimaries(prim);

  EPartInfo result{0.,0.,0.,0.,0,0};
  result.Ehit = Ehit;
  result.Esamp = Esamp;

  for (size_t i=0; i<prim.size(); i++) {
    const G4ParticleDefinition* pd = prim[i]->GetParticleDefinition();
    G4double wt = prim[i]->GetWeight();
    G4double E = prim[i]->GetKineticEnergy();

    if (G4CMP::IsPhonon(pd)) {
      result.Nphon += wt;
      result.Ephon += E*wt;
    }

    if (G4CMP::IsChargeCarrier(pd)) {
      result.Nchg += wt;
      result.Echg += E*wt;
    }
  }

  std::for_each(prim.begin(), prim.end(), Delete<G4PrimaryParticle>());
  prim.clear();

  return result;
}


// Main test is here

int main(int argc, char* argv[]) {
  if (argc < 4) {
    G4cerr << "Usage: " << argv[0] << " <Ehit> <Esamp> <Lattice>" << G4endl
	   << "\tEnergies should be in eV" << G4endl;
    ::exit(1);
  }

  G4double Ehit = strtod(argv[1],NULL) * eV;
  G4double Esamp = strtod(argv[2],NULL) * eV;
  G4String lname = argv[3];
  G4String mname = "G4_"+lname;

  G4int verbose = (argc>4) ? atoi(argv[4]) : 0;

  // MUST USE 'new', SO THAT G4SolidStore CAN DELETE
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4Tubs* crystal = new G4Tubs("GeCrystal", 0., 5.*cm, 1.*cm, 0., 360.*deg);
  G4LogicalVolume* lv = new G4LogicalVolume(crystal, mat, crystal->GetName());
  G4PVPlacement* pv = new G4PVPlacement(0, G4ThreeVector(), lv, lv->GetName(),
					0, false, 1);

  lattice = G4LatticeManager::Instance()->LoadLattice(pv,lname);

  partition = new G4CMPEnergyPartition(mat, lattice);
  partition->SetVerboseLevel(verbose);

  G4double bandgap = lattice->GetBandGapEnergy();

  G4double E=0., Esum=0., Esum2=0.;	// Sum and sum of squares for RMS
  G4double N, Nsum=0., Nsum2=0.;

  G4int nTest=1000;
  for (G4int i=0; i<nTest; i++) {
    EPartInfo iTest = testPartition(Ehit, Esamp);

    if (verbose) {
      G4cout << i << " : " << iTest.Nchg << " e/h " << iTest.Echg/keV << " keV"
	     << " + bandgaps " << iTest.Nchg*bandgap/2./keV << " keV"
	     << " " << iTest.Nphon << " phonons " << iTest.Ephon/keV << " keV"
	     << G4endl;
    }

    E = iTest.Echg+iTest.Ephon; Esum += E; Esum2 += E*E;
    N = iTest.Nchg+iTest.Nphon; Nsum += N; Nsum2 += N*N;
  }

  Esum /= nTest; Esum2 /= nTest;
  Nsum /= nTest; Nsum2 /= nTest;

  G4double E_SD = sqrt(fabs(Esum2 - Esum*Esum));
  G4double N_SD = sqrt(fabs(Nsum2 - Nsum*Nsum));

  G4cout << " Ehit " << Ehit/keV << " keV "
	 << "   <N> = " << Nsum << " +/- " << N_SD
	 << "   <E> = " << Esum/keV << " +/- " << E_SD << " keV" << G4endl;
}
