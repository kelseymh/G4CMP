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

#include "G4LatticeReader.hh"
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
  : verboseLevel(vb), psLatfile(0), pLattice(0), fMapPath(""),
    fToken(""), fValue(0.), fMap(""), fsPol(""), fPol(-1), fNX(0), fNY(0),
    f3Vec(0.,0.,0.), fLastBasis(-1),
    fDataDir(getenv("G4LATTICEDATA") ?
      (const char*)getenv("G4LATTICEDATA"):"./CrystalMaps"),
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
  fLastBasis = -1;			// Reset index counter for basis vectors

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

  pLattice->SetBasis();	// Fill or complete right-handed basis vectors
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

  // Extract path from filename to use in finding .ssv map files
  size_t lastdir = filepath.last('/');
  if (lastdir == std::string::npos) fMapPath = ".";	// No path at all
  else fMapPath = filepath(0,lastdir);

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
  if (fToken == "vdir")     return ProcessNMap();	// Direction vector map
  if (fToken == "vg")       return ProcessMap();	// Velocity magnitudes
  if (fToken == "dyn")      return ProcessConstants();	// Dynamical parameters
  if (fToken == "emass")    return ProcessMassTensor();	// e- mass eigenvalues
  if (fToken == "basis")    return ProcessBasisVector(); // Crystal axis dir
  if (fToken == "valley")   return ProcessEulerAngles(fToken); // e- drift dirs
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

// Read unit vector components for basis vectors (in order, b1, b2, b3)

G4bool G4LatticeReader::ProcessBasisVector() {
  *psLatfile >> f3Vec;
  if (verboseLevel>1) G4cout << " ProcessBasisVector " << f3Vec << G4endl;

  ++fLastBasis;
  if (fLastBasis>2) {
    G4cerr << " ERROR too many basis vectors.  Ignorning " << fLastBasis
	   << G4endl;
    return false;
  }

  pLattice->SetBasis(fLastBasis, f3Vec.unit());
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

// Read map filename, polarization, and binning dimensions 

G4bool G4LatticeReader::ReadMapInfo() {
  *psLatfile >> fMap >> fsPol >> fNX >> fNY;
  if (verboseLevel>1)
    G4cout << " ReadMapInfo " << fMap << " " << fsPol
	   << " " << fNX << " " << fNY << G4endl;

  if (fNX < 0 || fNX >= G4LatticeLogical::MAXRES) {
    G4cerr << "G4LatticeReader: Invalid map theta dimension " << fNX << G4endl;
    return false;
  }

  if (fNY < 0 || fNY >= G4LatticeLogical::MAXRES) {
    G4cerr << "G4LatticeReader: Invalid map phi dimension " << fNY << G4endl;
    return false;
  }

  // Prepend path to data files to map filename
  fMap = fMapPath + "/" + fMap;

  // Convert string code (L,ST,LT) to polarization index
  fsPol.toLower();
  fPol = ( (fsPol=="l")  ? 0 :		// Longitudinal
	   (fsPol=="st") ? 1 :		// Slow-transverse
	   (fsPol=="ft") ? 2 :		// Fast-transverse
	   -1 );			// Invalid code

  if (fPol<0 || fPol>2) {
    G4cerr << "G4LatticeReader: Invalid polarization code " << fsPol << G4endl;
    return false;
  }

  return true;
}

G4bool G4LatticeReader::ProcessMap() {
  if (!ReadMapInfo()) {		// Get specific parameters for map to load
    G4cerr << "G4LatticeReader: Unable to process mapfile directive." << G4endl;
    return false;
  }

  return pLattice->LoadMap(fNX, fNY, fPol, fMap);
}

G4bool G4LatticeReader::ProcessNMap() {
  if (!ReadMapInfo()) {		// Get specific parameters for map to load
    G4cerr << "G4LatticeReader: Unable to process mapfile directive." << G4endl;
    return false;
  }

  return pLattice->Load_NMap(fNX, fNY, fPol, fMap);
}
