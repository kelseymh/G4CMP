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
//
// 20210818  Expand test results to cover Fano disabled, downsampling, add
//		output reporting full range (min to max) of track counts.
// 20210820  Add bias voltage by hand to test Luke energy estimator.
// 20220914  G4CMP-322 -- Include N_SD in reporting output.

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

EPartInfo testPartition(G4double Ehit) {
  partition->DoPartition(Ehit, 0.);

  std::vector<G4PrimaryParticle*> prim;
  partition->GetPrimaries(prim);

  EPartInfo result{Ehit, g4cmp->GetSamplingEnergy(), 0., 0., 0, 0};

  for (size_t i=0; i<prim.size(); i++) {
    const G4ParticleDefinition* pd = prim[i]->GetParticleDefinition();
    G4double wt = prim[i]->GetWeight();
    G4double E = prim[i]->GetKineticEnergy();

    if (G4CMP::IsPhonon(pd)) {
      result.Nphon++;
      result.Ephon += E*wt;
    }

    if (G4CMP::IsChargeCarrier(pd)) {
      result.Nchg++;
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

  g4cmp->SetSamplingEnergy(Esamp);
  g4cmp->SetLukeSampling(-1.);			// Let partitioner do scaling

  partition = new G4CMPEnergyPartition(pv);
  partition->SetVerboseLevel(verbose);
  partition->SetBiasVoltage(50.*volt);		// For testing Luke sampling

  G4double bandgap = lattice->GetBandGapEnergy();

  G4double E=0., Esum=0., Esum2=0.;	// Sum and sum of squares for RMS
  G4double N, Nsum=0., Nsum2=0., Nchg=0., Nchg2=0.;
  G4double Nmin=DBL_MAX, Nmax=-DBL_MAX;

  G4int nTest=1000;
  for (G4int i=0; i<nTest; i++) {
    EPartInfo iTest = testPartition(Ehit);

    if (verbose) {
      G4cout << i << "\t" << iTest.Nchg << " e/h " << iTest.Echg/keV << " keV"
	     << " + bandgaps " << iTest.Nchg*bandgap/2./keV << " keV"
	     << " " << iTest.Nphon << " phonons " << iTest.Ephon/keV << " keV"
	     << G4endl;
    }

    E = iTest.Echg+iTest.Ephon; Esum += E; Esum2 += E*E;
    N = iTest.Nchg+iTest.Nphon; Nsum += N; Nsum2 += N*N;
    Nchg += iTest.Nchg; Nchg2 += iTest.Nchg*iTest.Nchg;

    if (iTest.Nchg<Nmin) Nmin = iTest.Nchg;
    if (iTest.Nchg>Nmax) Nmax = iTest.Nchg;
  }

  Esum /= nTest; Esum2 /= nTest;
  Nsum /= nTest; Nsum2 /= nTest;
  Nchg /= nTest; Nchg2 /= nTest;

  G4double E_SD = sqrt(fabs(Esum2 - Esum*Esum));
  G4double N_SD = sqrt(fabs(Nsum2 - Nsum*Nsum));
  G4double Q_SD = sqrt(fabs(Nchg2 - Nchg*Nchg));

  G4double Nexp = 2.*Ehit/lattice->GetPairProductionEnergy();
  G4double Nsig = (g4cmp->FanoStatisticsEnabled() 
		   ? sqrt(Nexp*lattice->GetFanoFactor()) : 0.);

  // With downsampling, the number expecting should be reduced by scale
  if (Esamp > 0 && Ehit > Esamp) {
    Nexp *= Esamp/Ehit;
    Nsig *= Esamp/Ehit;
  }

  G4cout << " Ehit " << Ehit/keV << " keV expect "
	 << Nexp << " +/- " << Nsig << " tracks (N_SD = " << N_SD << ")"
	 << " (" << lattice->GetLattice()->GetName()
	 << " Epair " << lattice->GetPairProductionEnergy()/eV << " eV)"
	 << "\n                 <Nchg> = " << Nchg << " +/- " << Q_SD
	 << " range " << Nmin << " to " << Nmax
	 << "\n <E> = " << Esum/keV << " +/- " << E_SD << " keV"
	 << " + bandgaps " << Nchg*bandgap/2./keV
	 << " = " << (Esum+Nchg*bandgap/2.)/keV << " keV " << G4endl;
}
