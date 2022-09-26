/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Usage: testTemperature <T> <Material>
//
// Specify temperature (K) and lattice directory/material name.
// Geant4 material will be set as "G4_<Material>".
//
// NOTE: Typical devices should be run at 0.1, 0.05, 0.02 K.
//
// 20220921  G4CMP-319 -- Exercise new temperature parameter and distribution.

#include "globals.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include <stdlib.h>
#include <string>
#include <vector>


// Global variables for use in tests

namespace {
  G4CMPConfigManager* g4cmp = G4CMPConfigManager::Instance();
  G4LatticePhysical* lattice = 0;
}


// Construct lattice with a "fake" physical volume

G4LatticePhysical* BuildLattice(const G4String& lname) {
  G4String mname = "G4_"+lname;

  // MUST USE 'new', SO THAT G4SolidStore CAN DELETE
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4Tubs* crystal = new G4Tubs(lname+"Crystal", 0., 5.*cm, 1.*cm, 0., 360.*deg);
  G4LogicalVolume* lv = new G4LogicalVolume(crystal, mat, crystal->GetName());
  G4PVPlacement* pv = new G4PVPlacement(0, G4ThreeVector(), lv, lv->GetName(),
					0, false, 1);

  lattice = G4LatticeManager::Instance()->LoadLattice(pv,lname);

  return lattice;
}


void DrawPDF(G4double temp, G4int verbose) {
  if (!verbose) return;

  G4double kT = k_Boltzmann*temp;	// Convert temperature to energy

  std::vector<std::string> graph(20, std::string(75,' '));
  for (G4int i=0; i<75; i++) {
    G4double Ei = kT * i/7.5;
    G4double pdf = G4CMP::MaxwellBoltzmannPDF(temp, Ei);
    G4double pdfj = std::min(int((1.-pdf)*20),19);	// Rows from top down
    if (verbose>1) {
      G4cout << " col " << i << " (" << Ei/eV << " eV) has "
	     << " pdf " << pdf << " = row " << pdfj << G4endl;
    }

    graph[pdfj][i] = '*';
  }
  
  G4cout << "\nMaxwell-Boltzmann PDF for " << temp/kelvin << " K" << G4endl;
  for (const auto& line : graph) G4cout << " | " << line << G4endl;
  G4cout << " +-" << std::string(75,'-') << G4endl;
}

void ThrowEnergies(G4double temp, G4int verbose) {
  G4double kT = k_Boltzmann*temp;	// Convert temperature to energy

  G4double Ebin = kT/7.5;		// 75 bins from 0 to 10*kT
  std::map<G4int,G4int> energy;		// Primitive histogram, energy bins

  G4int Nthrow = 50000;
  for (G4int i=0; i<Nthrow; i++) {	// Throw energy points to fill hist
    G4double E = G4CMP::ChooseThermalEnergy(temp);
    G4double j = int(E/Ebin);
    energy[j]++;			// Creates new bin, increments
    if (verbose>1) {
      G4cout << " Throw " << i << " E = " << E/eV << " eV,"
	     << " bin " << j << " has " << energy[j] << " entries"
	     << G4endl;
    }
  }

  // Generate text-based histogram of what got chosen
  G4int iMax = -1;
  for (const auto& Ei: energy) iMax = std::max(iMax, Ei.second);

  G4int yScale = iMax/19;		// Rescale so max bin fills plot

  std::vector<std::string> graph(20, std::string(75,' '));
  for (const auto& Ei: energy) {
    G4int i = Ei.first; if (i>74) continue;
    G4int n = Ei.second;
    G4int j = n/yScale;
    G4double jFrac = double(n%yScale)/yScale;

    if (verbose>1) {
      G4cout << " col " << i << " (" << (i+0.5)*Ebin/eV << " eV) has "
	     << " N " << n << " = row " << j << G4endl;
    }

    graph[19-j][i] = (jFrac<0.5 ? '.' : ':');
  }
  
  G4cout << "\n" << Nthrow << " random energies for " << temp/kelvin << " K"
	 << G4endl;
  for (const auto& line : graph) G4cout << " | " << line << G4endl;
  G4cout << " +-" << std::string(75,'-') << G4endl;

}

G4double IntegratePDF(G4double temp) {
  G4double kT = k_Boltzmann*temp;	// Convert temperature to energy

  G4double sum = 0., dE = kT/20.;
  for (G4double Ei=0.; Ei<20.*kT; Ei+=dE) {
    sum += G4CMP::MaxwellBoltzmannPDF(temp, Ei+dE/2.) / 20.;
  }

  return sum;
}

// Main test is here

int main(int argc, char* argv[]) {
  if (argc < 3) {
    G4cerr << "Usage: " << argv[0] << " <Temperature> <Material>" << G4endl
	   << "\tTemperature should be in kelvins" << G4endl;
    ::exit(1);
  }

  G4double temp = strtod(argv[1],NULL) * kelvin;
  G4String lname = argv[2];
  G4int verbose = (argc>3) ? atoi(argv[3]) : 0;

  BuildLattice(lname);
  lattice->SetTemperature(temp);
  if (verbose) lattice->Dump(G4cout);

  G4double kT = k_Boltzmann*temp;	// Convert temperature to energy

  // Draw a text picture of the PDF from 0 to 10 kT, 20 rows x 60 values
  DrawPDF(temp, verbose);
  ThrowEnergies(temp, verbose);

  // Calculate a few points, make sure lattice temperature is set
  G4double pdfTemp = G4CMP::MaxwellBoltzmannPDF(temp, kT);
  G4double pdfLat = G4CMP::MaxwellBoltzmannPDF(lattice->GetTemperature(), kT);
  G4double pdfHalf = G4CMP::MaxwellBoltzmannPDF(temp, kT/2.);
  G4double pdfDouble = G4CMP::MaxwellBoltzmannPDF(temp, kT*2.);
  G4double pdfFive = G4CMP::MaxwellBoltzmannPDF(temp, kT*5.);
  G4double pdfTen = G4CMP::MaxwellBoltzmannPDF(temp, kT*10.);

  G4cout << "T = " << temp/kelvin << " K,"
	 << " E = kT = " << kT/eV << " eV" << G4endl
	 << "PDF(E)   " << pdfTemp << " from lattice " << pdfLat << G4endl
	 << "PDF(E/2) " << pdfHalf << " PDF(2E) " << pdfDouble << G4endl
	 << "PDF(5E)  " << pdfFive << " PDF(10E) " << pdfTen << G4endl;

  // Estimate integral of PDF from 0. up to 20*kT, should be unity
  G4double sum = IntegratePDF(temp);
  G4cout << "Integral [0,20kT] = " << sum << G4endl
	 << G4endl;

  // Result is success if equal PDFs and unit integral
  return (pdfTemp == pdfLat && fabs(sum-1.) < 1e-6);
}
