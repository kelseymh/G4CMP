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
    mElectron(electron_mass_c2/c_squared) {}

G4LatticeReader::~G4LatticeReader() {
  delete psLatfile; psLatfile = 0;
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
       if (name == "beta")   pLattice->SetBeta(fValue);
  else if (name == "gamma")  pLattice->SetGamma(fValue);
  else if (name == "lambda") pLattice->SetLambda(fValue);
  else if (name == "mu")     pLattice->SetMu(fValue);
  else if (name == "scat")   pLattice->SetScatteringConstant(fValue*s*s*s);
  else if (name == "b")      pLattice->SetScatteringConstant(fValue*s*s*s);
  else if (name == "decay")  pLattice->SetAnhDecConstant(fValue*s*s*s*s);
  else if (name == "a")      pLattice->SetAnhDecConstant(fValue*s*s*s*s);
  else if (name == "ldos")   pLattice->SetLDOS(fValue);
  else if (name == "stdos")  pLattice->SetSTDOS(fValue);
  else if (name == "ftdos")  pLattice->SetFTDOS(fValue);
  else if (name == "vsound") pLattice->SetSoundSpeed(fValue*m/s);
  else if (name == "escat")  pLattice->SetElectronScatter(fValue*m);
  else if (name == "l0_e")   pLattice->SetElectronScatter(fValue*m);
  else if (name == "hscat")  pLattice->SetHoleScatter(fValue*m);
  else if (name == "l0_h")   pLattice->SetHoleScatter(fValue*m);
  else if (name == "hmass")  pLattice->SetHoleMass(fValue*mElectron);
  else if (name == "ivfield") pLattice->SetIVField(fValue);	// in V/m
  else if (name == "ivrate") pLattice->SetIVRate(fValue/s);
  else if (name == "ivpower") pLattice->SetIVExponent(fValue);
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
  if (verboseLevel>1)
    G4cout << " ProcessConstants " << beta << " " << gamma
	   << " " << lambda << " " << mu << G4endl;

  pLattice->SetDynamicalConstants(beta, gamma, lambda, mu);
  return psLatfile->good();
}


// Read lattice constants and angles for specified symmetry

G4bool G4LatticeReader::ProcessCrystalGroup(const G4String& name) {
  // Input buffers for reading; different crystals need different data
  G4double a=0., b=0., c=0., alpha=0., beta=0., gamma=0.;
  G4String unit;

  G4CMPCrystalGroup::Bravais group = G4CMPCrystalGroup::Group(name);
  switch (group) {
  case G4CMPCrystalGroup::amorphous:
    a=b=c=1./angstrom; break;			// No lattice constants
  case G4CMPCrystalGroup::cubic:
    *psLatfile >> a; b=c=a; break;		// Equal sides, orthogonal
  case G4CMPCrystalGroup::tetragonal:
  case G4CMPCrystalGroup::hexagonal:
    *psLatfile >> a >> c; b=c; break;		// Two sides, orthogonal
  case G4CMPCrystalGroup::orthorhombic:
    *psLatfile >> a >> b >> c; break;		// Three sides, orthogonal
  case G4CMPCrystalGroup::rhombohedral:
    *psLatfile >> a >> alpha >> unit; b=c=a; break;
  case G4CMPCrystalGroup::monoclinic:
    *psLatfile >> a >> b >> c >> alpha >> unit; break;
  case G4CMPCrystalGroup::triclinic:
    *psLatfile >> a >> b >> c >> alpha >> beta >> gamma >> unit; break;
  default: break;
  }

  G4double degOrRad = unit.empty() ? 0. : G4UnitDefinition::GetValueOf(unit);
  pLattice->SetCrystal(group, a*angstrom, b*angstrom, c*angstrom,
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

  pLattice->SetCij(p-1,q-1,value);	// Convention is C11-C66,
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
  G4String unit;
  *psLatfile >> phi >> theta >> psi >> unit;
  if (verboseLevel>1)
    G4cout << " ProcessEulerAngles " << name << " " << phi << " " 
	   << theta << " " << psi << " " << unit << G4endl;

  if (name != "valley") {
    G4cerr << "G4LatticeReader: Unknown rotation matrix " << name << G4endl;
    return false;
  }

  G4double degOrRad = G4UnitDefinition::GetValueOf(unit);
  pLattice->AddValley(phi*degOrRad, theta*degOrRad, psi*degOrRad);
  return psLatfile->good();
}
