// testNRyield: Exercise different Lindhard-yield (NIEL) functions
//
// Usage: testNRyield Emin Emax unit [Material = G4_Si]
//
// Arguments: Emin, Emax specify energy range to test, with the specified
//            units (e.g., "1 10 keV").  Optional Material should be a
//            Geant4 NIST name ("G4_something"), used to get both the
//            projectile and target (Z,A) values.
//
// Code will generate 50 steps for the given energy range, and call each of
// the available NIEL functions (this utility must be updated when new NIEL
// functions are added to the library).  The ionization yield from each
// function will be printed in a tab-delimited format suitable for plotting.
//
// 20250102  Michael Kelsey

#include "globals.hh"
#include "G4CMPLewinSmithNIEL.hh"
#include "G4CMPLindhardNIEL.hh"
#include "G4CMPImpactTunlNIEL.hh"
#include "G4CMPSarkisNIEL.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4VNIELPartition.hh"
#include <cmath>
#include <float.h>
#include <stdlib.h>


void testNRyield(G4double Emin, G4double Emax, const G4String& unit,
		 const G4String& material="G4_Si") {
  // Make sure that units and material strings are valid
  G4bool goodInput = true;
  if (!(goodInput &= G4UnitDefinition::IsUnitDefined(unit)))
    G4cerr << "ERROR: Invalid units '" << unit << "'" << G4endl;

  if (!(goodInput &= (G4UnitDefinition::GetCategory(unit) == "Energy")))
    G4cerr << "ERROR: Energy units (e.g., MeV) must be specified." << G4endl;

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* target = nist->FindOrBuildMaterial(material,true,true);
  goodInput &= (target != 0);

  // From material, Get Z and A for use as projectile (internal nuclear recoil)
  G4double Zin = target->GetZ(), Ain = target->GetA();

  // Sanity check: compare above with G4VNIELPartition functions
  // ==> Why are we getting errors from G4CMPImpactTunlNIEL for G4_Si???
  G4double Zniel = G4VNIELPartition::GetEffectiveZ(target);
  G4double Aniel = G4VNIELPartition::GetEffectiveA(target);

  if (!(goodInput &= (fabs(Zin-Zniel) < 0.5))) {
    G4cerr << "ERROR: Material::GetZ " << Zin << " doesn't match Zniel "
	   << Zniel << G4endl;
  }

  if (!(goodInput &= (fabs(Ain-Aniel) < 0.5))) {
    G4cerr << "ERROR: Material::GetA " << Ain << " doesn't match Aniel "
	   << Aniel << G4endl;
  }

  if (!goodInput) ::exit(1);		// If anything failed, abort

  // Instantiate single instance of each of the named yield functions
  const char* useNIEL[] = { "Lindhard", "LewinSmith", "Sarkis", "Impact" };
  const size_t nNIEL = sizeof(useNIEL)/sizeof(char*);

  // NOTE: Can't use G4CMPConfigManager to do this mapping, because it
  //       deletes the previous pointer when a new one is requested.
  const G4VNIELPartition* NIELfunc[nNIEL] = {
    new G4CMPLindhardNIEL, new G4CMPLewinSmithNIEL,
    new G4CMPSarkisNIEL, new G4CMPImpactTunlNIEL };

  // Output will be tab-delimited columns for all the NIEL functions
  G4cout << "Energy";
  for (size_t i=0; i<nNIEL; i++) G4cout << "\t" << useNIEL[i];
  G4cout << G4endl;

  // Compute uniform logarithmic steps to see full range of yield values
  const G4int nStep = 50;
  const G4double Estep = (Emax-Emin)/nStep;

  // Loop over energy steps, get yield values at each energy
  G4double unitVal = G4UnitDefinition::GetValueOf(unit);

  for (G4int iStep=0; iStep<=nStep; iStep++) {	// Avoids cumulative rounding
    G4double E = Emin + iStep*Estep;
    G4cout << E;

    for (size_t iNIEL=0; iNIEL<nNIEL; iNIEL++) {
      G4double Y =  NIELfunc[iNIEL]->PartitionNIEL(E*unitVal, target, Zin, Ain);
      G4cout << "\t" << Y;
    }
    G4cout << G4endl;
  }
}


// MAIN PROGRAM

int main(int argc, char* argv[]) {
  // Get required arguments, or report usage and exit
  if (argc < 4) {
    G4cerr << "Usage: testNRyield Emin Emax unit [Material = G4_Si]\n\n"
	   << "Arguments: Emin, Emax: specify energy range to test\n"
	   << "           units: Valid Geant4 unit string (e.g., 'MeV')\n"
	   << "           material: Valid Geant4 NIST name; default is G4_Si\n"
	   << "\nCode generates tab-delimited table of steps in energy,\n"
	   << "computing ionization yield Y from each of the NIEL-partition\n"
	   << "(Lindhard model) functions available in G4CMP.  The material\n"
	   << "is used for both target and projectile, corresponding to\n"
	   << "internal nuclear recoils in a detector."
	   << G4endl;
    ::exit(1);
  }

  G4double Emin = strtod(argv[1],0);
  G4double Emax = strtod(argv[2],0);
  G4String unit = argv[3];
  G4String mat  = (argc>4) ? argv[4] : "G4_Si";

  G4cout << "Nuclear recoils in " << mat << " " << Emin << "-" << Emax
	 << " " << unit << G4endl;

  testNRyield(Emin, Emax, unit, mat);
}
