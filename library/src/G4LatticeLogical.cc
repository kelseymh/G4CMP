/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/src/G4LatticeLogical.cc
/// \brief Implementation of the G4LatticeLogical class
//
// $Id$
//
// 20140218  Add new charge-carrier parameters to output
// 20140306  Allow valley filling using Euler angles directly
// 20140318  Compute electron mass scalar (Herring-Vogt) from tensor
// 20140324  Include inverse mass-ratio tensor
// 20140408  Add valley momentum calculations
// 20140425  Add "effective mass" calculation for electrons
// 20150601  Add mapping from electron velocity back to momentum
// 20160517  Add basis vectors for lattice, to use with Miller orientation
// 20160520  Add reporting function to format valley Euler angles
// 20160614  Add elasticity tensors and density (set from G4Material) 
// 20160624  Add direct calculation of phonon kinematics from elasticity
// 20160627  Interpolate values from lookup tables
// 20160629  Add post-constuction initialization (for tables, computed pars)
// 20160630  Drop loading of K-Vg lookup table files
// 20160701  Add interface to set elements of reduced elasticity matrix
// 20160727  Store Debye energy for phonon primaries, support different access
// 20170523  Add interface for axis vector of valleys
// 20170525  Add "rule of five" copy/move semantics
// 20170527  Drop unnecessary <fstream>
// 20170810  Add parameters for IV scattering matrix terms
// 20170821  Add support for separate D0 and D1 optical deformation potentials

#include "G4LatticeLogical.hh"
#include "G4CMPPhononKinematics.hh"	// **** THIS BREAKS G4 PORTING ****
#include "G4CMPPhononKinTable.hh"	// **** THIS BREAKS G4 PORTING ****
#include "G4CMPConfigManager.hh"	// **** THIS BREAKS G4 PORTING ****
#include "G4CMPUnitsTable.hh"		// **** THIS BREAKS G4 PORTING ****
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeLogical::G4LatticeLogical(const G4String& name)
  : verboseLevel(0), fName(name), fDensity(0.), fNImpurity(0.),
    fPermittivity(1.),
    fElasticity{}, fElReduced{}, fHasElasticity(false),
    fpPhononKin(0), fpPhononTable(0),
    fA(0), fB(0), fLDOS(0), fSTDOS(0), fFTDOS(0),
    fBeta(0), fGamma(0), fLambda(0), fMu(0),
    fVSound(0.), fVTrans(0.), fL0_e(0.), fL0_h(0.), 
    mElectron(electron_mass_c2/c_squared),
    fHoleMass(mElectron), fElectronMass(mElectron), fElectronMDOS(mElectron),
    fBandGap(0.), fPairEnergy(0.), fFanoFactor(1.),
    fMassTensor(G4Rep3x3(mElectron,0.,0.,0.,mElectron,0.,0.,0.,mElectron)),
    fMassInverse(G4Rep3x3(1/mElectron,0.,0.,0.,1/mElectron,0.,0.,0.,1/mElectron)),
    fAlpha(0.), fAcDeform(0.), fIVField(0.), fIVRate(0.), fIVExponent(0.) {
  for (G4int i=0; i<G4PhononPolarization::NUM_MODES; i++) {
    for (G4int j=0; j<KVBINS; j++) {
      for (G4int k=0; k<KVBINS; k++) {
	fKVMap[i][j][k].set(0.,0.,0.);
      }
    }
  }
}

G4LatticeLogical::~G4LatticeLogical() {
  delete fpPhononKin; fpPhononKin = 0;
  delete fpPhononTable; fpPhononTable = 0;
}

// Copy and move operators (to handle owned pointers)

G4LatticeLogical::G4LatticeLogical(const G4LatticeLogical& rhs)
  : G4LatticeLogical() { *this = rhs; }

G4LatticeLogical::G4LatticeLogical(G4LatticeLogical&& rhs)
  : G4LatticeLogical() { std::swap(*this, rhs); }

G4LatticeLogical& G4LatticeLogical::operator=(const G4LatticeLogical& rhs) {
  if (this == &rhs) return *this;	// Avoid unnecessary work;

  verboseLevel = rhs.verboseLevel;
  fName = rhs.fName;
  fCrystal = rhs.fCrystal;
  std::copy(rhs.fBasis, rhs.fBasis+3, fBasis);
  fDensity = rhs.fDensity;
  fNImpurity = rhs.fNImpurity;
  fPermittivity = rhs.fPermittivity;
  fHasElasticity = rhs.fHasElasticity;
  fA = rhs.fA;
  fB = rhs.fB;
  fLDOS = rhs.fLDOS;
  fSTDOS = rhs.fSTDOS;
  fFTDOS = rhs.fFTDOS;
  fBeta = rhs.fBeta;
  fGamma = rhs.fGamma;
  fLambda = rhs.fLambda;
  fMu = rhs.fMu;
  fDebye = rhs.fDebye;
  fVSound = rhs.fVSound;
  fVTrans = rhs.fVTrans;
  fL0_e = rhs.fL0_e;
  fL0_h = rhs.fL0_h;
  fHoleMass = rhs.fHoleMass;
  fElectronMass = rhs.fElectronMass;
  fElectronMDOS = rhs.fElectronMDOS;
  fBandGap = rhs.fBandGap;
  fPairEnergy = rhs.fPairEnergy;
  fFanoFactor = rhs.fFanoFactor;
  fMassTensor = rhs.fMassTensor;
  fMassInverse = rhs.fMassInverse;
  fMassRatioSqrt = rhs.fMassRatioSqrt;
  fMInvRatioSqrt = rhs.fMInvRatioSqrt;
  fValley = rhs.fValley;
  fValleyAxis = rhs.fValleyAxis;
  fAlpha = rhs.fAlpha;
  fAcDeform = rhs.fAcDeform;
  fIVDeform = rhs.fIVDeform;
  fIVEnergy = rhs.fIVEnergy;
  fIVField = rhs.fIVField;
  fIVRate = rhs.fIVRate;
  fIVExponent = rhs.fIVExponent;

  if (!rhs.fpPhononKin)   fpPhononKin = new G4CMPPhononKinematics(this);
  if (!rhs.fpPhononTable) fpPhononTable = new G4CMPPhononKinTable(fpPhononKin);

  SetElReduced(rhs.fElReduced);
  FillElasticity();

  for (G4int i=0; i<G4PhononPolarization::NUM_MODES; i++) {
    for (G4int j=0; j<KVBINS; j++) {
      for (G4int k=0; k<KVBINS; k++) {
	fKVMap[i][j][k] = rhs.fKVMap[i][j][k];
      }
    }
  }

  return *this;
}

G4LatticeLogical& G4LatticeLogical::operator=(G4LatticeLogical&& rhs) {
  std::swap(*this, rhs);
  return *this;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
//Copy or set components of reduced elasticity matrix
/////////////////////////////////////////////////////////////
void G4LatticeLogical::SetElReduced(const ReducedElasticity& mat) {
  for (size_t i=0; i<6; i++) {
    for (size_t j=0; j<6; j++) {
      fElReduced[i][j] = mat[i][j];
    }
  }

  fHasElasticity = true;
}

void G4LatticeLogical::SetCpq(G4int p, G4int q, G4double value) {
  if (p>0 && p<7 && q>0 && q<7) fElReduced[p-1][q-1] = value;
  fHasElasticity = true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
// Configure crystal symmetry group and lattice spacing/angles
/////////////////////////////////////////////////////////////
void G4LatticeLogical::SetCrystal(G4CMPCrystalGroup::Bravais group, G4double a,
				  G4double b, G4double c, G4double alpha,
				  G4double beta, G4double gamma) {
  fCrystal.Set(group, alpha, beta, gamma);	// Defines unit cell axes

  fBasis[0] = a*fCrystal.axis[0];	// Basis vectors include spacing
  fBasis[1] = b*fCrystal.axis[1];
  fBasis[2] = c*fCrystal.axis[2];
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
//Configured derived parameters and tables after loading
/////////////////////////////////////////////////////////////
void G4LatticeLogical::Initialize(const G4String& newName) {
  if (!newName.empty()) SetName(newName);

  CheckBasis();				// Ensure complete, right handed frame

  // If elasticity matrix available, create phonon calculator
  if (fHasElasticity) {
    FillElasticity();			// Unpack reduced matrix to full Cijkl
    if (!fpPhononKin) fpPhononKin = new G4CMPPhononKinematics(this);
  }

  /***** USE OUR OWN INTERPOLATION, THIS IS TOO SLOW
  if (fpPhononKin) fpPhononTable = new G4CMPPhononKinTable(fpPhononKin);
  *****/

  // Populate phonon lookup tables if not read from files
  FillMaps();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
//Complete basis vectors: right-handed, possibly orthonormal
/////////////////////////////////////////////////////////////
void G4LatticeLogical::CheckBasis() {
  static const G4ThreeVector origin(0.,0.,0.);
  if (fBasis[0].isNear(origin,1e-6)) fBasis[0].set(1.,0.,0.);
  if (fBasis[1].isNear(origin,1e-6)) fBasis[1].set(0.,1.,0.);
  if (fBasis[2].isNear(origin,1e-6)) fBasis[2] = fBasis[0].cross(fBasis[1]);

  // Ensure that all basis vectors are unit
  fBasis[0].setMag(1.);
  fBasis[1].setMag(1.);
  fBasis[2].setMag(1.);

  if (fBasis[0].cross(fBasis[1]).dot(fBasis[2]) < 0.) {
    G4cerr << "ERROR G4LatticeLogical has a left-handed basis!" << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Unpack reduced elasticity tensor into full four-dimensional Cijkl
void G4LatticeLogical::FillElasticity() {
  fCrystal.FillElReduced(fElReduced);		// Apply symmetry conditions
  /* The reduced matrix looks like this:
   * Cxxxx, Cxxyy, Cxxzz, Cxxyz, Cxxxz, Cxxxy
   * Cxxyy, Cyyyy, Cyyzz, Cyyyz, Cyyxz, Cyyxy
   * Cxxyy, Cyyzz, Czzzz, Czzyz, Czzxz, Czzxy
   * Cxxyz, Cyyyz, Czzyz, Cyzyz, Cyzxz, Cyzxy
   * Cxxxz, Cyyxz, Czzxz, Cyzxz, Cxzxz, Cxzxy
   * Cxxxy, Cyyxy, Czzxy, Cyzxy, Cxzxy, Cxyxy
   */

  auto reducedIdx = [](size_t i, size_t j) -> size_t {
    // i == j -> i
    // i == 0 && j == 1 -> 5
    // i == 0 && j == 2 -> 4
    // i == 1 && j == 2 -> 3
    if (i > 2 || j > 2) {
      G4Exception("G4LatticeLogical::FillElasticity",
                  "Lattice011",
                  EventMustBeAborted,
                  "Indices can only span 0 to 2 (x to z).");
    }

    if (i == j) return i;
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 5;
    if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
    if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;

    // Should never get to this point (Exception thrown above
    return 0;
  };

  /* This could potentially be sped up because of various symmetries:
   * Cijkl = Cjikl = Cijlk
   * Cxxyy = Cyyxx
   * But it probably doesn't matter because this is only done at init anyway.
   * Plus it may get slowed down by messing with the CPU's branch prediction
   * and compiler loop unrolling.
   */
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          fElasticity[i][j][k][l] = fElReduced[reducedIdx(i, j)][reducedIdx(k, l)];
        }
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Populate lookup tables using kinematics calculator

void G4LatticeLogical::FillMaps() {
  if (!fpPhononKin) return;			// Can't fill without solver

  G4ThreeVector k;
  for (G4int itheta = 0; itheta<KVBINS; itheta++) {
    G4double theta = itheta*pi/(KVBINS-1);	// Last entry is at pi

    for (G4int iphi = 0; iphi<KVBINS; iphi++) {
      G4double phi = iphi*twopi/(KVBINS-1);	// Last entry is at 2pi

      k.setRThetaPhi(1.,theta,phi);
      for (G4int mode=0; mode<G4PhononPolarization::NUM_MODES; mode++) {
	fKVMap[mode][itheta][iphi] = fpPhononKin->getGroupVelocity(mode,k);
      }
    }
  }

  if (verboseLevel) {
    G4cout << "G4LatticeLogical::FillMaps populated " << KVBINS
	   << " bins in theta and phi for all polarizations." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon group velocity vector

G4ThreeVector G4LatticeLogical::MapKtoVg(G4int polarizationState,
					 const G4ThreeVector& k) const {
  return ( (fpPhononKin && G4CMPConfigManager::UseKVSolver())
	   ? ComputeKtoVg(polarizationState,k)
	   : LookupKtoVg(polarizationState,k) );
}

G4ThreeVector G4LatticeLogical::ComputeKtoVg(G4int polarizationState,
					     const G4ThreeVector& k) const {  
  if (!fpPhononKin) {
    G4Exception("G4LatticeLogical::ComputeKtoV", "Lattice001",
		RunMustBeAborted, "Phonon kinematics not available.");
  }

  return fpPhononKin->getGroupVelocity(polarizationState,k);
}

G4ThreeVector G4LatticeLogical::LookupKtoVg(G4int polarizationState,
					    const G4ThreeVector& k) const {  
  if (fpPhononTable)
    return (fpPhononTable->interpGroupVelocity(polarizationState, k.unit())*
	    fpPhononTable->interpGroupVelocity_N(polarizationState, k.unit()).unit()
	    );

  G4int iTheta, iPhi;		// Bin indices
  G4double dTheta, dPhi;	// Offsets in bin for interpolation
  if (!FindLookupBins(k, iTheta, iPhi, dTheta, dPhi)) {
    G4Exception("G4LatticeLogical::LookupKtoVDir", "Lattice006",
		EventMustBeAborted, "Interpolation failed.");
    return G4ThreeVector();
  }

  /**** Returns direct bin value
  const G4ThreeVector& vdir = fKVMap[polarizationState][iTheta][iPhi];
  ****/

  // Bilinear interpolation using the four corner bins (i,j) to (i+1,j+1)
  G4ThreeVector vdir =
    ( (1.-dTheta)*(1.-dPhi)*fKVMap[polarizationState][iTheta][iPhi] +
      dTheta*(1.-dPhi)*fKVMap[polarizationState][iTheta+1][iPhi] +
      (1.-dTheta)*dPhi*fKVMap[polarizationState][iTheta][iPhi+1] +
      dTheta*dPhi*fKVMap[polarizationState][iTheta+1][iPhi+1] );

  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::MapKtoVDir theta,phi="
	   << k.theta() << " " << k.phi()
	   << " : ith,iph " << iTheta << " " << iPhi
	   << " : dir " << vdir << G4endl;
  }

  return vdir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Get theta, phi bins and offsets for interpolation

G4bool 
G4LatticeLogical::FindLookupBins(const G4ThreeVector& k,
				 G4int& iTheta, G4int& iPhi,
				 G4double& dTheta, G4double& dPhi) const {
  G4double tStep = pi/(KVBINS-1);	// Last element is upper edge
  G4double pStep = twopi/(KVBINS-1);

  G4double theta = k.getTheta();	// Normalize theta to [0,pi)
  if (theta<0) theta+=pi;

  G4double phi = k.getPhi();		// Normalize phi to [0,twopi)
  if (phi<0) phi += twopi;

  dTheta = theta/tStep;
  iTheta = int(dTheta);
  dTheta -= iTheta;			// Fraction of bin width

  dPhi = phi/pStep;
  iPhi = int(dPhi);
  dPhi -= iPhi;				// Fraction of bin width

  return (iTheta<KVBINS && iPhi<KVBINS);	// Sanity check on bin indexing
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Convert electron momentum to valley velocity, wavevector, and HV vector

G4ThreeVector 
G4LatticeLogical::MapPtoV_el(G4int ivalley, const G4ThreeVector& p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoV_el " << ivalley << " " << p_e
	   << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return vToN.inverse()*(GetMInvTensor()*(vToN*p_e/c_light));
}

G4ThreeVector 
G4LatticeLogical::MapV_elToP(G4int ivalley, const G4ThreeVector& v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToP " << ivalley << " " << v_e
	   << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return vToN.inverse()*(GetMassTensor()*(vToN*v_e*c_light));
}

G4ThreeVector
G4LatticeLogical::MapV_elToK_HV(G4int ivalley, const G4ThreeVector &v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToK_HV " << ivalley << " " << v_e
     << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return GetSqrtInvTensor()*GetMassTensor()*vToN*v_e/hbar_Planck;
}

G4ThreeVector 
G4LatticeLogical::MapPtoK_valley(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoK " << ivalley << " " << p_e
	   << G4endl;

  p_e /= hbarc;					// Convert to wavevector
  return p_e.transform(GetValley(ivalley));	// Rotate into valley frame
}

G4ThreeVector 
G4LatticeLogical::MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoK_HV " << ivalley << " " << p_e
	   << G4endl;

  p_e.transform(GetValley(ivalley));		// Rotate into valley frame
  return GetSqrtInvTensor() * p_e/hbarc;	// Herring-Vogt transformation
}

G4ThreeVector 
G4LatticeLogical::MapK_HVtoK_valley(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoK_valley " << ivalley << " " << k_HV
	   << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley
  return k_HV;
}

G4ThreeVector
G4LatticeLogical::MapK_HVtoK(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoK " << ivalley << " " << k_HV
     << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley
  k_HV.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  return k_HV;
}

G4ThreeVector 
G4LatticeLogical::MapK_HVtoP(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoP " << ivalley << " " << k_HV
	   << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley 
  k_HV.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  k_HV *= hbarc;			// Convert wavevector to momentum
  return k_HV;
}

G4ThreeVector 
G4LatticeLogical::MapK_valleyToP(G4int ivalley, G4ThreeVector k) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_valleyToP " << ivalley << " " << k
	   << G4endl;

  k.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  k *= hbarc;				// Convert wavevector to momentum
  return k;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Apply energy-momentum relationship for electron transport

G4double  
G4LatticeLogical::MapPtoEkin(G4int iv, G4ThreeVector p) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoEkin " << iv << " " << p << G4endl;

  p.transform(GetValley(iv));			// Rotate to valley frame

  // Compute kinetic energy component by component, then sum
  return (0.5/c_squared) * (p.x()*p.x()*fMassInverse.xx() +
			    p.y()*p.y()*fMassInverse.yy() +
			    p.z()*p.z()*fMassInverse.zz());
}

G4double
G4LatticeLogical::MapV_elToEkin(G4int iv, G4ThreeVector v) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToEkin " << iv << " " << v << G4endl;

  v.transform(GetValley(iv));			// Rotate to valley frame

  // Compute kinetic energy component by component, then sum
  return 0.5 * (v.x()*v.x()*fMassTensor.xx() +
          v.y()*v.y()*fMassTensor.yy() +
          v.z()*v.z()*fMassTensor.zz());
}

// Compute effective "scalar" electron mass to match energy/momentum relation

G4double 
G4LatticeLogical::GetElectronEffectiveMass(G4int iv,
					   const G4ThreeVector& p) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::GetElectronEffectiveMass " << iv
	   << " " << p << G4endl;

  return 0.5*p.mag2()/c_squared/MapPtoEkin(iv,p);	// Non-relativistic
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Store electron mass tensor using diagonal elements

void G4LatticeLogical::SetMassTensor(G4double mXX, G4double mYY, G4double mZZ) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::SetMassTensor " << mXX << " " << mYY
	   << " " << mZZ << " *m_e" << G4endl;
  }

  // NOTE:  Use of G4RotationMatrix not appropriate here, as matrix is
  //        not normalized.  But CLHEP/Matrix not available in GEANT4.
  fMassTensor.set(G4Rep3x3(mXX*mElectron, 0., 0.,
			   0., mYY*mElectron, 0.,
			   0., 0., mZZ*mElectron));

  FillMassInfo();
}

void G4LatticeLogical::SetMassTensor(const G4RotationMatrix& etens) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::SetMassTensor " << etens << G4endl;
  }

  // Check if mass tensor already has electron mass, or is just coefficients
  G4bool hasEmass = (etens.xx()/mElectron > 1e-3 ||
		     etens.yy()/mElectron > 1e-3 ||
		     etens.zz()/mElectron > 1e-3);
  G4double mscale = hasEmass ? 1. : mElectron;

  // NOTE:  Use of G4RotationMatrix not appropriate here, as matrix is
  //        not normalized.  But CLHEP/Matrix not available in GEANT4.
  fMassTensor.set(G4Rep3x3(etens.xx()*mscale, 0., 0.,
			   0., etens.yy()*mscale, 0.,
			   0., 0., etens.zz()*mscale));

  FillMassInfo();
}

// Compute derived quantities from user-input mass tensor

void G4LatticeLogical::FillMassInfo() {
  // Herring-Vogt scalar mass of electron, used for sound speed calcs  
  fElectronMass = 3. / ( 1./fMassTensor.xx() + 1./fMassTensor.yy()
			 + 1./fMassTensor.zz() );  

  // Density of states effective mass, used for intervalley scattering
  fElectronMDOS = cbrt(fMassTensor.xx()*fMassTensor.yy()*fMassTensor.zz());

  // 1/m mass tensor used for k and v calculations in valley coordinates
  fMassInverse.set(G4Rep3x3(1./fMassTensor.xx(), 0., 0.,
			    0., 1./fMassTensor.yy(), 0.,
			    0., 0., 1./fMassTensor.zz()));

  // Mass ratio tensor used for scattering and field calculations
  fMassRatioSqrt.set(G4Rep3x3(sqrt(fMassTensor.xx()/fElectronMass), 0., 0.,
			      0., sqrt(fMassTensor.yy()/fElectronMass), 0.,
			      0., 0., sqrt(fMassTensor.zz()/fElectronMass)));

  fMInvRatioSqrt.set(G4Rep3x3(1./fMassRatioSqrt.xx(), 0., 0.,
			      0., 1./fMassRatioSqrt.yy(), 0.,
			      0., 0., 1./fMassRatioSqrt.zz()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Store drifting-electron valley using Euler angles

void G4LatticeLogical::AddValley(G4double phi, G4double theta, G4double psi) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::AddValley " << psi << " " << theta
	   << " " << psi << " rad" << G4endl;
  }

  // Extend vector first, then fill last value, to reduce temporaries
  fValley.resize(fValley.size()+1);
  fValley.back().set(phi,theta,psi);

  // NOTE:  Rotation matrices take external vector along valley axis to X-hat
  fValleyAxis.push_back(fValley.back().inverse()*G4ThreeVector(1.,0.,0.));
}

// Store rotation matrix and corresponding axis vector for valley

void G4LatticeLogical::AddValley(const G4RotationMatrix& valley) {
  fValley.push_back(valley);

  // NOTE:  Rotation matrices take external vector along valley axis to X-hat
  fValleyAxis.push_back(valley.inverse()*G4ThreeVector(1.,0.,0.));
}

// Transform for drifting-electron valleys in momentum space

const G4RotationMatrix& G4LatticeLogical::GetValley(G4int iv) const {
  if (verboseLevel>1) G4cout << "G4LatticeLogical::GetValley " << iv << G4endl;

  if (iv >=0 && iv < (G4int)NumberOfValleys()) return fValley[iv];

  if (verboseLevel)
    G4cerr << "G4LatticeLogical ERROR: No such valley " << iv << G4endl;
  return G4RotationMatrix::IDENTITY;
}

const G4ThreeVector& G4LatticeLogical::GetValleyAxis(G4int iv) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::GetValleyAxis " << iv << G4endl;

  if (iv >=0 && iv < (G4int)NumberOfValleys()) return fValleyAxis[iv];

  if (verboseLevel)
    G4cerr << "G4LatticeLogical ERROR: No such valley " << iv << G4endl;

  static const G4ThreeVector nullVec(0.,0.,0.);
  return nullVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Set Debye energy for phonon partitioning from alternative parameters

void G4LatticeLogical::SetDebyeFreq(G4double nu) { fDebye = nu*h_Planck; }

void G4LatticeLogical::SetDebyeTemp(G4double temp) { fDebye = temp*k_Boltzmann;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump structure in format compatible with reading back

void G4LatticeLogical::Dump(std::ostream& os) const {
  os << "# " << fName << " crystal lattice parameters"
     << "\n# density " << fDensity/(g/cm3) << " g/cm3"
     << std::endl;

  if (fHasElasticity) DumpCrystalInfo(os);

  os << "# Phonon propagation parameters"
     << "\ndyn " << fBeta/GPa << " " << fGamma/GPa  << " "
     << fLambda/GPa  << " " << fMu/GPa << " GPa"
     << "\nscat " << fB/s3 << " s3" << " decay " << fA/s4 << " s4" 
     << "\nLDOS " << fLDOS << " STDOS " << fSTDOS << " FTDOS " << fFTDOS
     << "\nDebye " << fDebye/eV << " eV"
     << std::endl;

  os << "# Charge carrier propagation parameters"
     << "\nbandgap " << fBandGap/eV << " eV"
     << "\npairEnergy " << fPairEnergy/eV << " eV"
     << "\nfanoFactor " << fFanoFactor
     << "\nvsound " << fVSound/(m/s) << " m/s"
     << "\nl0_e " << fL0_e/um << " um"
     << "\nl0_h " << fL0_h/um << " um"
     << std::endl;

  os << "# Charge carrier masses [m(electron) units]"
     << "\nhmass " << fHoleMass/mElectron
     << "\nemass " << fMassTensor.xx()/mElectron
     << " " << fMassTensor.yy()/mElectron
     << " " << fMassTensor.zz()/mElectron << std::endl;

  os << "# Inverse mass tensor: " << fMassInverse.xx()*mElectron
     << " " << fMassInverse.yy()*mElectron
     << " " << fMassInverse.zz()*mElectron
     << " * 1/m(electron)" << std::endl
     << "# Herring-Vogt scalar mass: " << fElectronMass/mElectron << std::endl
     << "# Density of states mass: " << fElectronMDOS/mElectron << std::endl
     << "# sqrt(tensor/scalar): " << fMassRatioSqrt.xx()
     << " " << fMassRatioSqrt.yy()
     << " " << fMassRatioSqrt.zz()
     << std::endl;

  for (size_t i=0; i<NumberOfValleys(); i++) {
    DumpValley(os, i);
  }

  os << "# Intervalley scattering parameters"
     << "\nalpha " << fAlpha*eV << " /eV"
     << "\nepsilon " << fPermittivity
     << "\nneutDens " << fNImpurity * cm3 << " /cm3"
     << "\nacDeform " << fAcDeform/eV << " eV"
     << "\n ivDeform "; DumpList(os, fIVDeform, "eV/cm");
  os << "\n ivEnergy "; DumpList(os, fIVEnergy, "eV");
  os << std::endl;

  os << "# Edelweiss intervalley scattering parameters"
     << "\nivField " << fIVField/(volt/m) << " V/m"
     << "\nivRate " << fIVRate/hertz << " Hz"
     << "\nivPower " << fIVExponent << std::endl;
}

// Print out Euler angles of requested valley

void G4LatticeLogical::DumpValley(std::ostream& os, G4int iv) const {
  if (iv < 0 || static_cast<size_t>(iv) >= NumberOfValleys()) return;

  os << "valley " << fValley[iv].phi()/deg
     << " " << fValley[iv].theta()/deg
     << " " << fValley[iv].psi()/deg
     << " deg" << std::endl;
}

// Print out crystal symmetry information

void G4LatticeLogical::DumpCrystalInfo(std::ostream& os) const {
  G4double a=fBasis[0].mag()/angstrom;		// Lattice params for printing
  G4double b=fBasis[1].mag()/angstrom;
  G4double c=fBasis[2].mag()/angstrom;

  os << fCrystal.Name() << " ";
  switch (fCrystal.group) {		// Lattice constants and angles
  case G4CMPCrystalGroup::cubic:
    os << a << " Ang"; break;
  case G4CMPCrystalGroup::tetragonal:
  case G4CMPCrystalGroup::hexagonal:
    os << a << c << " Ang"; break;
  case G4CMPCrystalGroup::orthorhombic:
    os << a << " " << b << " " << c << " Ang"; break;
  case G4CMPCrystalGroup::rhombohedral:
    os << a << " Ang " << fCrystal.alpha()/deg << " deg"; break;
  case G4CMPCrystalGroup::monoclinic:
    os << a << " " << b << " " << c << " Ang "
       << fCrystal.alpha()/deg << " deg"; break;
  case G4CMPCrystalGroup::triclinic:
    os << a << " " << b << " " << c << " Ang "
       << fCrystal.alpha()/deg << " " << fCrystal.beta()/deg << " "
       << fCrystal.gamma()/deg << " deg"; break;
  default: break;
  }
  os << std::endl;

  switch (fCrystal.group) {		// Reduced elasticity tensor
  case G4CMPCrystalGroup::tetragonal:
    DumpCpq(os,1,6);				// Plus all below
  case G4CMPCrystalGroup::hexagonal:
    DumpCpq(os,1,3); DumpCpq(os,3,3); DumpCpq(os,6,6);	// Plus all below
  case G4CMPCrystalGroup::cubic:
    DumpCpq(os,4,4);				// Plus all below
  case G4CMPCrystalGroup::amorphous:
    DumpCpq(os,1,1); DumpCpq(os,1,2); break;
  case G4CMPCrystalGroup::rhombohedral:
    DumpCpq(os,1,1); DumpCpq(os,1,2); DumpCpq(os,1,3);
    DumpCpq(os,1,4); DumpCpq(os,1,5);
    DumpCpq(os,3,3); DumpCpq(os,4,4); DumpCpq(os,6,6); break;
  case G4CMPCrystalGroup::monoclinic:
    DumpCpq(os,1,6); DumpCpq(os,2,6); DumpCpq(os,3,6);	// Plus all below
    DumpCpq(os,4,5);
  case G4CMPCrystalGroup::orthorhombic:
    for (int p=1; p<7; p++) {		// Upper corner and lower diagonal
      for (int q=p; q<4; q++) DumpCpq(os,p,q);
      if (p>3) DumpCpq(os,p,p);
    }
    break;
  case G4CMPCrystalGroup::triclinic:	// Entire upper half, no zeroes
    for (int p=1; p<7; p++) for (int q=p; q<7; q++) DumpCpq(os,p,q);
    break;
  default: break;
  }
}

// Print out elasticity tensor element with units

void G4LatticeLogical::DumpCpq(std::ostream& os, G4int p, G4int q) const {
  os << "Cpq " << p << " " << q << " " << GetCpq(p,q)/GPa << " GPa"
     << std::endl;
}

// Print out list of values scaled by specified unit

void G4LatticeLogical::DumpList(std::ostream& os,
				const std::vector<G4double>& vlist,
				const G4String& unit) const {
  if (vlist.empty()) return;		// Avoid unnecessary work

  G4double uval = G4UnitDefinition::GetValueOf(unit);
  for (size_t i=0; i<vlist.size(); i++) {
    os << vlist[i]/uval << " ";
  }

  os << unit;
}
