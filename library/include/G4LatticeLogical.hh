/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/include/G4LatticeLogical.hh
/// \brief Definition of the G4LatticeLogical class
//
// $Id$
//
// 20131114  Add verbosity for diagnostic output
// 20131115  Expose maximum array dimensions for use by LatticeReader
// 20140218  Add support for charge-carrier functionality
// 20140306  Allow valley filling using Euler angles directly
// 20140313  Allow electron mass filling with diagonal elements
// 20140319  Add "extra" mass tensors with precomputed expressions
// 20140324  Add intervalley scattering parameters
// 20140408  Add valley momentum calculations
// 20140425  Add "effective mass" calculation for electrons
// 20150601  Add mapping from electron velocity back to momentum
// 20160517  Add basis vectors for lattice, to use with Miller orientation
// 20160520  Add reporting function to format valley Euler angles
// 20160614  Add elasticity tensors and density (set from G4Material) 
// 20160624  Add direct calculation of phonon kinematics from elasticity
// 20160629  Add post-constuction initialization (for tables, computed pars)
// 20160630  Drop loading of K-Vg lookup table files

#ifndef G4LatticeLogical_h
#define G4LatticeLogical_h

#include "globals.hh"
#include "G4CMPCrystalGroup.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PhononPolarization.hh"
#include <iosfwd>
#include <vector>

class G4CMPPhononKinematics;
class G4CMPPhononKinTable;

// Arrays for full and reduced elasticity matrices
class G4LatticeLogical {
public:
  typedef G4double Elasticity[3][3][3][3];
  typedef G4double ReducedElasticity[6][6];

public:
  G4LatticeLogical(const G4String& name="");
  virtual ~G4LatticeLogical();

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  void SetName(const G4String& name) { fName = name; }
  const G4String& GetName() const { return fName; }

  // Compute derived quantities, fill tables, etc. after setting parameters
  void Initialize(const G4String& name="");

  // Dump structure in format compatible with reading back
  void Dump(std::ostream& os) const;

  // Get group velocity magnitude, direction for input polarization and wavevector
  // NOTE:  Wavevector must be in lattice symmetry frame (X == symmetry axis)
  virtual G4ThreeVector MapKtoVg(G4int pol, const G4ThreeVector& k) const;

  virtual G4double MapKtoV(G4int pol, const G4ThreeVector& k) const {
    return MapKtoVg(pol,k).mag();
  }

  virtual G4ThreeVector MapKtoVDir(G4int pol, const G4ThreeVector& k) const {
    return MapKtoVg(pol,k).unit();
  }

  // Convert between electron momentum and valley velocity or HV wavevector
  // NOTE:  Input vector must be in lattice symmetry frame (X == symmetry axis)
  // NOTE:  Pass by value below to avoid creating temporary vectors
  G4ThreeVector MapPtoV_el(G4int ivalley, const G4ThreeVector& p_e) const;
  G4ThreeVector MapV_elToP(G4int ivalley, const G4ThreeVector& v_el) const;
  G4ThreeVector MapV_elToK_HV(G4int ivalley, const G4ThreeVector& v_el) const;
  G4ThreeVector MapPtoK_valley(G4int ivalley, G4ThreeVector p_e) const;
  G4ThreeVector MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const;
  G4ThreeVector MapK_HVtoP(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_HVtoK_valley(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_HVtoK(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_valleyToP(G4int ivalley, G4ThreeVector k) const;

  // Apply energy relationships for electron transport
  G4double MapPtoEkin(G4int ivalley, G4ThreeVector p_e) const;
  G4double MapV_elToEkin(G4int ivalley, G4ThreeVector v_e) const;

  // Configure crystal symmetry group and lattice spacing/angles
  void SetCrystal(G4CMPCrystalGroup::Bravais group, G4double a, G4double b,
		  G4double c, G4double alpha, G4double beta, G4double gamma);

  // Get specified basis vector (returns null if invalid index)
  const G4ThreeVector& GetBasis(G4int i) const {
    static const G4ThreeVector nullVec(0.,0.,0.);
    return (i>=0 && i<3 ? fBasis[i] : nullVec);
  }

  // Physical parameters of lattice (density, elasticity)
  void SetDensity(G4double val) { fDensity = val; }

  G4double GetDensity() const { return fDensity; }
  const Elasticity& GetElasticity() const { return fElasticity; }
  G4double GetCijkl(G4int i, G4int j, G4int k, G4int l) const {
    return fElasticity[i][j][k][l];
  }

  void SetElReduced(const ReducedElasticity& mat);
  const ReducedElasticity& GetElReduced() const { return fElReduced; }

  void SetCij(G4int i, G4int j, G4double value);
  G4double GetCij(G4int i, G4int j) const { return fElReduced[i][j]; }

  // Parameters for phonon production and propagation
  void SetDynamicalConstants(G4double Beta, G4double Gamma,
			     G4double Lambda, G4double Mu) {
    fBeta=Beta; fGamma=Gamma; fLambda=Lambda; fMu=Mu;
  }

  void SetBeta(G4double Beta) { fBeta = Beta; }
  void SetGamma(G4double Gamma) { fGamma = Gamma; }
  void SetLambda(G4double Lambda) { fLambda = Lambda; }
  void SetMu(G4double Mu) { fMu = Mu; }
  void SetScatteringConstant(G4double b) { fB=b; }
  void SetAnhDecConstant(G4double a) { fA=a; }
  void SetLDOS(G4double LDOS) { fLDOS=LDOS; }
  void SetSTDOS(G4double STDOS) { fSTDOS=STDOS; }
  void SetFTDOS(G4double FTDOS) { fFTDOS=FTDOS; }
  void SetDebyeFreq(G4double nu) { fDebNu=nu; }

  G4double GetBeta() const { return fBeta; }
  G4double GetGamma() const { return fGamma; }
  G4double GetLambda() const { return fLambda; }
  G4double GetMu() const { return fMu; }
  G4double GetScatteringConstant() const { return fB; }
  G4double GetAnhDecConstant() const { return fA; }
  G4double GetLDOS() const { return fLDOS; }
  G4double GetSTDOS() const { return fSTDOS; }
  G4double GetFTDOS() const { return fFTDOS; }
  G4double GetDebyeFreq() const { return fDebNu; }

  // Parameters and structures for charge carrier transport
  void SetBandGapEnergy(G4double bg) { fBandGap = bg; }
  void SetPairProductionEnergy(G4double pp) { fPairEnergy = pp; }
  void SetFanoFactor(G4double f) { fFanoFactor = f; }
  void SetSoundSpeed(G4double v) { fVSound = v; }
  void SetHoleScatter(G4double l0) { fL0_h = l0; }
  void SetHoleMass(G4double hmass) { fHoleMass = hmass; }
  void SetElectronScatter(G4double l0) { fL0_e = l0; }
  void SetMassTensor(const G4RotationMatrix& etens);
  void SetMassTensor(G4double mXX, G4double mYY, G4double mZZ);

  G4double GetBandGapEnergy() const             { return fBandGap; }
  G4double GetPairProductionEnergy() const      { return fPairEnergy; }
  G4double GetFanoFactor() const                { return fFanoFactor; }
  G4double GetSoundSpeed() const                { return fVSound; }
  G4double GetHoleScatter() const               { return fL0_h; }
  G4double GetHoleMass() const                  { return fHoleMass; }
  G4double GetElectronScatter() const           { return fL0_e; }
  G4double GetElectronMass() const 		{ return fElectronMass; }
  const G4RotationMatrix& GetMassTensor() const { return fMassTensor; }
  const G4RotationMatrix& GetMInvTensor() const { return fMassInverse; }
  const G4RotationMatrix& GetSqrtTensor() const { return fMassRatioSqrt; }
  const G4RotationMatrix& GetSqrtInvTensor() const { return fMInvRatioSqrt; }

  // Compute "effective mass" for electron to preserve E/p relationship
  G4double GetElectronEffectiveMass(G4int iv, const G4ThreeVector& p) const;

  // Transform for drifting-electron valleys in momentum space
  void AddValley(const G4RotationMatrix& valley) { fValley.push_back(valley); }
  void AddValley(G4double phi, G4double theta, G4double psi);
  void ClearValleys() { fValley.clear(); }

  size_t NumberOfValleys() const { return fValley.size(); }
  const G4RotationMatrix& GetValley(G4int iv) const;

  // Print out Euler angles of requested valley
  void DumpValley(std::ostream& os, G4int iv) const;

  // Parameters for electron intervalley scattering
  void SetIVField(G4double v)    { fIVField = v; }
  void SetIVRate(G4double v)     { fIVRate = v; }
  void SetIVExponent(G4double v) { fIVExponent = v; }

  G4double GetIVField() const    { return fIVField; }
  G4double GetIVRate() const     { return fIVRate; }
  G4double GetIVExponent() const { return fIVExponent; }

private:
  void CheckBasis();	// Initialize or complete (via cross) basis vectors
  void FillElasticity();	// Unpack reduced Cij into full Cijlk
  void FillMaps();	// Populate lookup tables using kinematics calculator
  void FillMassInfo();	// Called from SetMassTensor() to compute derived forms

  // Get theta, phi bins and offsets for interpolation
  G4bool FindLookupBins(const G4ThreeVector& k, G4int& iTheta, G4int& iPhi,
			G4double& dTheta, G4double& dPhi) const;

  // Use lookup table to get group velocity for phonons
  G4ThreeVector LookupKtoVg(G4int pol, const G4ThreeVector& k) const;

  // Use direct calculation to get group velocity for phonons
  G4ThreeVector ComputeKtoVg(G4int pol, const G4ThreeVector& k) const;

private:
  G4int verboseLevel;			    // Enable diagnostic output
  G4String fName;			    // Name of lattice for messages

  G4CMPCrystalGroup fCrystal;		    // Symmetry group, axis unit vectors
  G4ThreeVector fBasis[3];		    // Basis vectors for Miller indices
  G4double fDensity;			    // Material density (natural units)
  Elasticity fElasticity;	    	    // Full 4D elasticity tensor
  ReducedElasticity fElReduced;		    // Reduced 2D elasticity tensor
  G4bool fHasElasticity;		    // Flag valid elasticity tensors
  G4CMPPhononKinematics* fpPhononKin;	    // Kinematics calculator with tensor
  G4CMPPhononKinTable* fpPhononTable;	    // Kinematics interpolator

  // map for group velocity vectors
  enum { KVBINS=315 };			    // K-Vg lookup table binning
  G4ThreeVector fKVMap[G4PhononPolarization::NUM_MODES][KVBINS][KVBINS];

  G4double fA;       // Scaling constant for Anh.Dec. mean free path
  G4double fB;       // Scaling constant for Iso.Scat. mean free path
  G4double fLDOS;    // Density of states for L-phonons
  G4double fSTDOS;   // Density of states for ST-phonons
  G4double fFTDOS;   // Density of states for FT-phonons
  G4double fBeta, fGamma, fLambda, fMu; // dynamical constants for material
  G4double fDebNu;  // Debye frequency

  G4double fVSound;	// Speed of sound (longitudinal phonon)
  G4double fL0_e;	// Scattering length for electrons
  G4double fL0_h;	// Scattering length for holes

  const G4double mElectron;	 // Free electron mass (without G4's c^2)
  G4double fHoleMass;		 // Effective mass of +ve carrier
  G4double fElectronMass;	 // Effective mass (scalar) of -ve carrier
  G4double fBandGap;	 // Minimum band gap energy
  G4double fPairEnergy;   // electron-hole pair production average energy
  G4double fFanoFactor;   // Fano factor (duh)
  G4RotationMatrix fMassTensor;	 // Full electron mass tensor
  G4RotationMatrix fMassInverse; // Inverse electron mass tensor (convenience)
  G4RotationMatrix fMassRatioSqrt;       // SQRT of tensor/scalar ratio
  G4RotationMatrix fMInvRatioSqrt;       // SQRT of scalar/tensor ratio
  std::vector<G4RotationMatrix> fValley; // Electron transport directions

  G4double fIVField;		 // Transverse field for intervalley scattering
  G4double fIVRate;		 // Scale factor for IV scattering MFP
  G4double fIVExponent;		 // Power law for E-field in IV scattering
};

// Write lattice structure to output stream

inline std::ostream& 
operator<<(std::ostream& os, const G4LatticeLogical& lattice) {
  lattice.Dump(os);
  return os;
}

#endif
