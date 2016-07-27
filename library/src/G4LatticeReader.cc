/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4LatticeReader.hh
/// \brief Implementation of the G4LatticeReader class
//
// NOTE:  This reader class for logical lattices should be moved to
//	  materials/ after the 10.0 release (and this comment removed).
// $Id$
//
// 20131106  M.Kelsey -- Add const to getenv() to avoid compiler warning.
// 20131112  Throw exception if input file fails.
// 20131115  Check file input arguments for maps for validity before use;
//		move ctor, dtor here; check stream pointer before closing.
// 20140218  Add support for charge-carrier functionality
// 20140306  Use Euler angles directly to fill electron valley
// 20140313  Use diagonal terms directly to fill electron mass tensor
// 20140314  Add charge-carrier propagation parameters
// 20140324  Add intervalley scattering parameters
// 20160517  Add basis vectors for lattice
// 20160615  Add elasticity tensor (cubic lattice only)
// 20160630  Drop loading of K-Vg lookup table files
// 20160701  Withdraw seting basis vectors, set crystal symmetry instead
// 20160727  Require units in config file for all parameters,
//		define additional units for solid state physics use

#include "G4CMPConfigManager.hh"
#include "G4LatticeReader.hh"
#include "G4CMPCrystalGroup.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeLogical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include <limits>
#include <stdlib.h>


// Constructor and destructor

G4LatticeReader::G4LatticeReader(G4int vb)
  : verboseLevel(vb), psLatfile(0), pLattice(0),
    fToken(""), fValue(0.), f3Vec(0.,0.,0.),
    fDataDir(G4CMPConfigManager::GetLatticeDir()),
    mElectron(electron_mass_c2/c_squared) {
  DefineUnits();
}

G4LatticeReader::~G4LatticeReader() {
  delete psLatfile; psLatfile = 0;
}


// Define special time units for scattering rate coefficients

void G4LatticeReader::DefineUnits() {
  // Phonon scattering coeffient B [s^3]
  new G4UnitDefinition(     "second^3","s3"   ,"Time cubed",s*s*s);
  new G4UnitDefinition("millisecond^3","ms3"  ,"Time cubed",ms*ms*ms);
  new G4UnitDefinition("microsecond^3","mus3" ,"Time cubed",
		       microsecond*microsecond*microsecond);
  new G4UnitDefinition( "nanosecond^3","ns3"  ,"Time cubed",ns*ns*ns);
  new G4UnitDefinition( "picosecond^3","ps3"  ,"Time cubed",
			picosecond*picosecond*picosecond);
  new G4UnitDefinition(     "hertz^-3","Hz-3" ,"Time cubed",
			    1./(hertz*hertz*hertz));
  new G4UnitDefinition( "kilohertz^-3","kHz-3","Time cubed",
			1./(kilohertz*kilohertz*kilohertz));
  new G4UnitDefinition( "megahertz^-3","MHz-3","Time cubed",
			1./(megahertz*megahertz*megahertz));

  // Phonon anharmonic decay coefficient A [s^4]
  new G4UnitDefinition(     "second^4","s4"   ,"Time fourth",s*s*s*s);
  new G4UnitDefinition("millisecond^4","ms4"  ,"Time fourth",ms*ms*ms*ms);
  new G4UnitDefinition("microsecond^4","mus4" ,"Time fourth",
		       microsecond*microsecond*microsecond*microsecond);
  new G4UnitDefinition( "nanosecond^4","ns4"  ,"Time fourth",ns*ns*ns*ns);
  new G4UnitDefinition( "picosecond^4","ps4"  ,"Time fourth",
			picosecond*picosecond*picosecond*picosecond);
  new G4UnitDefinition(     "hertz^-4","Hz-4" ,"Time fourth",
			    1./(hertz*hertz*hertz*hertz));
  new G4UnitDefinition( "kilohertz^-4","kHz-4","Time fourth",
			1./(kilohertz*kilohertz*kilohertz*kilohertz));
  new G4UnitDefinition( "megahertz^-4","MHz-4","Time fourth",
			1./(megahertz*megahertz*megahertz*megahertz));

  // Stiffness (pressure) and frequency units suitable for solid state physics
  new G4UnitDefinition("gigapascal", "GPa", "Pressure",  1e9*pascal);
  new G4UnitDefinition( "terahertz", "THz", "Frequency", 1e12*hertz);

  // Velocity (?!? Why aren't these already defined in Geant4 ?!?)
  new G4UnitDefinition(     "meters/second",  "m/s", "Velocity", m/s);
  new G4UnitDefinition( "kilometers/second", "km/s", "Velocity", km/s);
  new G4UnitDefinition("millimeters/second", "mm/s", "Velocity", mm/s);
  new G4UnitDefinition("centimeters/second", "cm/s", "Velocity", cm/s);
}


// Main drivers to read configuration from file or stream

G4LatticeLogical* G4LatticeReader::MakeLattice(const G4String& filename) {
  if (verboseLevel) G4cout << "G4LatticeReader " << filename << G4endl;

  if (!OpenFile(filename)) {
    G4ExceptionDescription msg;
    msg << "Unable to open " << filename;
    G4Exception("G4LatticeReader::MakeLattice", "Lattice001",
		FatalException, msg);
    return 0;
  }

  pLattice = new G4LatticeLogical;	// Create lattice to be filled

  G4bool goodLattice = true;
  while (!psLatfile->eof()) {
    goodLattice &= ProcessToken();
  }
  CloseFile();

  if (!goodLattice) {
    G4ExceptionDescription msg;
    msg << "Error reading lattice from " << filename;
    G4Exception("G4LatticeReader::MakeLattice", "Lattice002",
		FatalException, msg);
    delete pLattice;
    pLattice = 0;
  }

  return pLattice;	// Lattice complete; return pointer with ownership
}


// Open local file or file found under data path

G4bool G4LatticeReader::OpenFile(const G4String& filename) {
  if (verboseLevel)
    G4cout << "G4LatticeReader::OpenFile " << filename << G4endl;

  G4String filepath = filename;
  psLatfile = new std::ifstream(filepath);
  if (!psLatfile->good()) {			// Local file not found
    filepath = fDataDir + "/" + filename;
    psLatfile->open(filepath);			// Try data directory
    if (!psLatfile->good()) {
      CloseFile();
      return false;
    }
    if (verboseLevel>1) G4cout << " Found file " << filepath << G4endl;
  }

  return true;
}

// Close and delete input stream

void G4LatticeReader::CloseFile() {
  if (psLatfile) psLatfile->close();
  delete psLatfile;
  psLatfile = 0;
}


// Read next token from file, use it to store next data into lattice

G4bool G4LatticeReader::ProcessToken() {
  fToken = "";
  *psLatfile >> fToken;
  if (fToken.empty() || psLatfile->eof()) return true;	// End of file reached

  if (verboseLevel>1) G4cout << " ProcessToken " << fToken << G4endl;

  fToken.toLower();
  if (fToken.contains('#')) return SkipComments();	// Ignore rest of line
  if (fToken == "dyn")      return ProcessConstants();	// Dynamical parameters
  if (fToken == "stiffness" ||
      fToken == "cij")      return ProcessStiffness();  // Elasticity element
  if (fToken == "emass")    return ProcessMassTensor();	// e- mass eigenvalues
  if (fToken == "valley")   return ProcessEulerAngles(fToken); // e- drift dirs
  if (G4CMPCrystalGroup::Group(fToken) >= 0)		// Crystal dimensions
                            return ProcessCrystalGroup(fToken);

  return ProcessValue(fToken);				// Single numeric value
}

// Eat remainder of line, assuming a '#' token was found

G4bool G4LatticeReader::SkipComments() {
  psLatfile->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  return true;		// Never fails
}

// Read double value from file, store based on name string

G4bool G4LatticeReader::ProcessValue(const G4String& name) {
  *psLatfile >> fValue;
  if (verboseLevel>1) G4cout << " ProcessValue " << fValue << G4endl;

  G4bool good = true;
       if (name == "beta")       pLattice->SetBeta(fValue*ProcessUnits("Pressure"));
  else if (name == "gamma")      pLattice->SetGamma(fValue*ProcessUnits("Pressure"));
  else if (name == "lambda")     pLattice->SetLambda(fValue*ProcessUnits("Pressure"));
  else if (name == "mu")         pLattice->SetMu(fValue*ProcessUnits("Pressure"));
  else if (name == "scat")       pLattice->SetScatteringConstant(fValue*ProcessUnits("Time cubed"));
  else if (name == "b")          pLattice->SetScatteringConstant(fValue*ProcessUnits("Time cubed"));
  else if (name == "decay")      pLattice->SetAnhDecConstant(fValue*ProcessUnits("Time fourth"));
  else if (name == "a")          pLattice->SetAnhDecConstant(fValue*ProcessUnits("Time fourth"));
  else if (name == "ldos")       pLattice->SetLDOS(fValue);
  else if (name == "stdos")      pLattice->SetSTDOS(fValue);
  else if (name == "ftdos")      pLattice->SetFTDOS(fValue);
  else if (name == "debyefreq")  pLattice->SetDebyeFreq(fValue*ProcessUnits("Frequency"));
  else if (name == "bandgap")    pLattice->SetBandGapEnergy(fValue*ProcessUnits("Energy"));
  else if (name == "pairenergy") pLattice->SetPairProductionEnergy(fValue*ProcessUnits("Energy"));
  else if (name == "fanofactor") pLattice->SetFanoFactor(fValue);
  else if (name == "vsound")     pLattice->SetSoundSpeed(fValue*ProcessUnits("Velocity"));
  else if (name == "escat")      pLattice->SetElectronScatter(fValue*ProcessUnits("Length"));
  else if (name == "l0_e")       pLattice->SetElectronScatter(fValue*ProcessUnits("Length"));
  else if (name == "hscat")      pLattice->SetHoleScatter(fValue*ProcessUnits("Length"));
  else if (name == "l0_h")       pLattice->SetHoleScatter(fValue*ProcessUnits("Length"));
  else if (name == "hmass")      pLattice->SetHoleMass(fValue*mElectron);
  else if (name == "ivfield")    pLattice->SetIVField(fValue*ProcessUnits("Electric field"));
  else if (name == "ivrate")     pLattice->SetIVRate(fValue*ProcessUnits("Frequency"));
  else if (name == "ivpower")    pLattice->SetIVExponent(fValue);
  else if (name == "ivexponent") pLattice->SetIVExponent(fValue);
  else {
    G4cerr << "G4LatticeReader: Unrecognized token " << name << G4endl;
    good = false;
  }

  return good;
}

G4bool G4LatticeReader::ProcessConstants() {
  G4double beta=0., gamma=0., lambda=0., mu=0.;
  *psLatfile >> beta >> gamma >> lambda >> mu;
  ProcessUnits("Pressure");

  if (verboseLevel>1)
    G4cout << " ProcessConstants " << beta << " " << gamma
	   << " " << lambda << " " << mu << " " << fUnitName << G4endl;

  pLattice->SetDynamicalConstants(beta*fUnits, gamma*fUnits, lambda*fUnits,
				  mu*fUnits);
  return psLatfile->good();
}


// Read lattice constants and angles for specified symmetry

G4bool G4LatticeReader::ProcessCrystalGroup(const G4String& name) {
  if (verboseLevel>1) G4cout << " ProcessCrystalGroup " << name << G4endl;


  // Input buffers for reading; different crystals need different data
  G4double a=0., b=0., c=0., alpha=0., beta=0., gamma=0.;
  G4double lunit=0., degOrRad=0.;		// Length and angle units

  G4CMPCrystalGroup::Bravais group = G4CMPCrystalGroup::Group(name);
  switch (group) {
  case G4CMPCrystalGroup::amorphous:
    a=b=c=1.; lunit=1.; break;			// No lattice constants
  case G4CMPCrystalGroup::cubic:
    *psLatfile >> a; b=c=a; 			// Equal sides, orthogonal
    lunit = ProcessUnits("Length");
    break;
  case G4CMPCrystalGroup::tetragonal:
  case G4CMPCrystalGroup::hexagonal:
    *psLatfile >> a >> c; b=c;			// Two sides, orthogonal
    lunit = ProcessUnits("Length");
    break;
  case G4CMPCrystalGroup::orthorhombic:
    *psLatfile >> a >> b >> c;			// Three sides, orthogonal
    lunit = ProcessUnits("Length");
    break;
  case G4CMPCrystalGroup::rhombohedral:
    *psLatfile >> a; b=c=a;
    lunit = ProcessUnits("Length");
    *psLatfile >> alpha;
    degOrRad = ProcessUnits("Angle");
    break;
  case G4CMPCrystalGroup::monoclinic:
    *psLatfile >> a >> b >> c;
    lunit = ProcessUnits("Length");
    *psLatfile >> alpha;
    degOrRad = ProcessUnits("Angle");
    break;
  case G4CMPCrystalGroup::triclinic:
    *psLatfile >> a >> b >> c;
    lunit = ProcessUnits("Length");
    *psLatfile >> alpha >> beta >> gamma;
    degOrRad = ProcessUnits("Angle");
    break;
  default: break;
  }

  if (verboseLevel>1) {
    G4cout << " a " << a*lunit/angstrom << " b " << b*lunit/angstrom
	   << " c " << c*lunit/angstrom << " Ang "
	   << " alpha " << alpha/deg << " beta " << beta/deg
	   << " gamma " << gamma/deg << " deg" << G4endl;
  }

  pLattice->SetCrystal(group, a*lunit, b*lunit, c*lunit,
		       alpha*degOrRad, beta*degOrRad, gamma*degOrRad);

  return psLatfile->good();
}

// Read element of reduced elasticity (stiffness) matrix

G4bool G4LatticeReader::ProcessStiffness() {
  G4int p=0, q=0;	// Indices of reduced matrix
  G4double value=0.;	// Matrix element in pascals

  *psLatfile >> p >> q >> value;
  if (verboseLevel>1)
    G4cout << "ProcessStiffness " << p << " " << q << " " << value << G4endl;

  // Convention for indices is C11-C66, but arguments below are (0,0) - (5,5)
  pLattice->SetCij(p-1,q-1,value*ProcessUnits("Pressure"));
  return psLatfile->good();
}

// Read diagonal scale factors for drift electron mass tensor

G4bool G4LatticeReader::ProcessMassTensor() {
  G4double mxx=1., myy=1., mzz=1.;
  *psLatfile >> mxx >> myy >> mzz;
  if (verboseLevel>1)
    G4cout << " ProcessMassTensor " << mxx << " " << myy << " " << mzz
	   << G4endl;

  pLattice->SetMassTensor(mxx, myy, mzz);
  return psLatfile->good();
}

// Read Euler angles (phi, theta, psi) for named rotation matrix

G4bool G4LatticeReader::ProcessEulerAngles(const G4String& name) {
  G4double phi=0., theta=0., psi=0.;
  *psLatfile >> phi >> theta >> psi;
  if (verboseLevel>1)
    G4cout << " ProcessEulerAngles " << name << " " << phi << " " 
	   << theta << " " << psi << G4endl;

  if (name != "valley") {
    G4cerr << "G4LatticeReader: Unknown rotation matrix " << name << G4endl;
    return false;
  }

  G4double degOrRad = ProcessUnits("Angle");
  pLattice->AddValley(phi*degOrRad, theta*degOrRad, psi*degOrRad);
  return psLatfile->good();
}


// Read expected dimensions for value from file, return scale factor
// Input argument "unitcat" may be comma-delimited list of categories

G4double G4LatticeReader::ProcessUnits(const G4String& unitcat) {
  *psLatfile >> fUnitName;
  return ProcessUnits(fUnitName, unitcat);
}

G4double G4LatticeReader::ProcessUnits(const G4String& unit,
				       const G4String& unitcat) {
  if (verboseLevel>1)
    G4cout << " ProcessUnits " << unit << " " << unitcat << G4endl;

  // Do processing -- invalid input string will cause fatal exception
  fUnitName = unit;
  fUnits    = G4UnitDefinition::GetValueOf(fUnitName);
  fUnitCat  = G4UnitDefinition::GetCategory(fUnitName);

  // If actual category doesn't match expected, throw exception
  if (!unitcat.contains(fUnitCat)) {
    G4ExceptionDescription msg;
    msg << "Expected " << unitcat << " units, got " << fUnitName << " ("
	<< fUnitCat << ")";
    G4Exception("G4LatticeReader::ProcessUnits", "Lattice003",
		FatalException, msg);
    return 0.;
  }

  return fUnits;	// Return value for convenient inlining
}
