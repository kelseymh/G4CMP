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
// 20160727  Store Debye energy for phonon primaries, support different access
// 20170523  Add interface for axis vector of valleys
// 20170525  Add "rule of five" copy/move semantics
// 20170810  Add parameters for IV scattering matrix terms
// 20170821  Add transverse sound speed, L->TT fraction
// 20170919  Add access to full lists of IV scattering matrix terms
// 20170928  Replace "pol" with "mode" for phonon states
// 20180815  F. Insulla -- Added IVRateQuad
// 20181001  M. Kelsey -- Clarify IV rate parameters systematically
// 20190704  M. Kelsey -- Add IV rate function selector for material
// 20190801  M. Kelsey -- Use G4ThreeVector buffer instead of pass-by-value,
//		precompute valley inverse transforms
// 20200608  Fix -Wshadow warnings from tempvec
// 20210919  M. Kelsey -- Allow SetVerboseLevel() from const instances.

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

  // Copy and move operators (to handle owned pointers)
  G4LatticeLogical(const G4LatticeLogical& rhs);
  G4LatticeLogical(G4LatticeLogical&& rhs);
  G4LatticeLogical& operator=(const G4LatticeLogical& rhs);
  G4LatticeLogical& operator=(G4LatticeLogical&& rhs);

  // Run-time configuration
  void SetVerboseLevel(G4int vb) const { verboseLevel = vb; }

  void SetName(const G4String& name) { fName = name; }
  const G4String& GetName() const { return fName; }

  // Compute derived quantities, fill tables, etc. after setting parameters
  void Initialize(const G4String& name="");

  // Dump structure in format compatible with reading back
  void Dump(std::ostream& os) const;

  // Get group velocity magnitude, direction for input polarization and wavevector
  // NOTE:  Wavevector must be in lattice symmetry frame (X == symmetry axis)
  virtual G4ThreeVector MapKtoVg(G4int mode, const G4ThreeVector& k) const;

  virtual G4double MapKtoV(G4int mode, const G4ThreeVector& k) const {
    return MapKtoVg(mode,k).mag();
  }

  virtual G4ThreeVector MapKtoVDir(G4int mode, const G4ThreeVector& k) const {
    return MapKtoVg(mode,k).unit();
  }

  // Convert between electron momentum and valley velocity or HV wavevector
  // NOTE:  Input vector must be in lattice symmetry frame (X == symmetry axis)
  G4ThreeVector MapPtoV_el(G4int ivalley, const G4ThreeVector& p_e) const;
  G4ThreeVector MapV_elToP(G4int ivalley, const G4ThreeVector& v_el) const;
  G4ThreeVector MapV_elToK_HV(G4int ivalley, const G4ThreeVector& v_el) const;
  G4ThreeVector MapPtoK_valley(G4int ivalley, const G4ThreeVector& p_e) const;
  G4ThreeVector MapPtoK_HV(G4int ivalley, const G4ThreeVector& p_e) const;
  G4ThreeVector MapK_HVtoP(G4int ivalley, const G4ThreeVector& k_HV) const;
  G4ThreeVector MapK_HVtoK_valley(G4int ivalley, const G4ThreeVector& k_HV) const;
  G4ThreeVector MapK_HVtoK(G4int ivalley, const G4ThreeVector& k_HV) const;
  G4ThreeVector MapK_valleyToP(G4int ivalley, const G4ThreeVector& k) const;

  // Apply energy relationships for electron transport
  G4double MapPtoEkin(G4int ivalley, const G4ThreeVector& p_e) const;
  G4double MapV_elToEkin(G4int ivalley, const G4ThreeVector& v_e) const;

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

  void SetImpurities(G4double val) { fNImpurity = val; }
  G4double GetImpurities() const { return fNImpurity; }

  void SetPermittivity(G4double val) { fPermittivity = val; }
  G4double GetPermittivity() const { return fPermittivity; }

  const Elasticity& GetElasticity() const { return fElasticity; }
  G4double GetCijkl(G4int i, G4int j, G4int k, G4int l) const {
    return fElasticity[i][j][k][l];
  }

  void SetElReduced(const ReducedElasticity& mat);
  const ReducedElasticity& GetElReduced() const { return fElReduced; }

  // Reduced elasticity tensor: C11-C66 interface for clarity
  void SetCpq(G4int p, G4int q, G4double value);
  G4double GetCpq(G4int p, G4int q) const { return fElReduced[p-1][q-1]; }

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
  void SetAnhTTFrac(G4double f) { fTTFrac=f; }
  void SetLDOS(G4double LDOS) { fLDOS=LDOS; }
  void SetSTDOS(G4double STDOS) { fSTDOS=STDOS; }
  void SetFTDOS(G4double FTDOS) { fFTDOS=FTDOS; }

  void SetDebyeEnergy(G4double energy) { fDebye = energy; }
  void SetDebyeFreq(G4double nu);
  void SetDebyeTemp(G4double temp);

  G4double GetBeta() const { return fBeta; }
  G4double GetGamma() const { return fGamma; }
  G4double GetLambda() const { return fLambda; }
  G4double GetMu() const { return fMu; }
  G4double GetScatteringConstant() const { return fB; }
  G4double GetAnhDecConstant() const { return fA; }
  G4double GetAnhTTFrac() const { return fTTFrac; }
  G4double GetLDOS() const { return fLDOS; }
  G4double GetSTDOS() const { return fSTDOS; }
  G4double GetFTDOS() const { return fFTDOS; }
  G4double GetDebyeEnergy() const { return fDebye; }

  // Parameters and structures for charge carrier transport
  void SetBandGapEnergy(G4double bg) { fBandGap = bg; }
  void SetPairProductionEnergy(G4double pp) { fPairEnergy = pp; }
  void SetFanoFactor(G4double f) { fFanoFactor = f; }
  void SetSoundSpeed(G4double v) { fVSound = v; }
  void SetTransverseSoundSpeed(G4double v) { fVTrans = v; }
  void SetHoleScatter(G4double l0) { fL0_h = l0; }
  void SetHoleMass(G4double hmass) { fHoleMass = hmass; }
  void SetElectronScatter(G4double l0) { fL0_e = l0; }
  void SetMassTensor(const G4RotationMatrix& etens);
  void SetMassTensor(G4double mXX, G4double mYY, G4double mZZ);

  G4double GetBandGapEnergy() const             { return fBandGap; }
  G4double GetPairProductionEnergy() const      { return fPairEnergy; }
  G4double GetFanoFactor() const                { return fFanoFactor; }
  G4double GetSoundSpeed() const                { return fVSound; }
  G4double GetTransverseSoundSpeed() const      { return fVTrans; }
  G4double GetHoleScatter() const               { return fL0_h; }
  G4double GetHoleMass() const                  { return fHoleMass; }
  G4double GetElectronScatter() const           { return fL0_e; }
  G4double GetElectronMass() const 		{ return fElectronMass; }
  G4double GetElectronDOSMass() const 		{ return fElectronMDOS; }
  const G4RotationMatrix& GetMassTensor() const { return fMassTensor; }
  const G4RotationMatrix& GetMInvTensor() const { return fMassInverse; }
  const G4RotationMatrix& GetSqrtTensor() const { return fMassRatioSqrt; }
  const G4RotationMatrix& GetSqrtInvTensor() const { return fMInvRatioSqrt; }

  // Compute "effective mass" for electron to preserve E/p relationship
  G4double GetElectronEffectiveMass(G4int iv, const G4ThreeVector& p) const;

  // Transform for drifting-electron valleys in momentum space
  void AddValley(const G4RotationMatrix& valley);
  void AddValley(G4double phi, G4double theta, G4double psi);
  void ClearValleys() {
    fValley.clear(); fValleyInv.clear();fValleyAxis.clear();
  }

  size_t NumberOfValleys() const { return fValley.size(); }
  const G4RotationMatrix& GetValley(G4int iv) const;
  const G4RotationMatrix& GetValleyInv(G4int iv) const;
  const G4ThreeVector& GetValleyAxis(G4int iv) const;

  // Print out Euler angles of requested valley
  void DumpValley(std::ostream& os, G4int iv) const;

  // Print out crystal symmetry information
  void DumpCrystalInfo(std::ostream& os) const;

  // Print out elasticity tensor element with units, in C11-C66 notation
  void DumpCpq(std::ostream& os, G4int p, G4int q) const;

  // Print out list of values, scaled by unit
  void DumpList(std::ostream& os, const std::vector<G4double>& vlist,
		const G4String& unit) const;

  // Parameters for electron intervalley scattering (Edelweiss, Linear, matrix)
  void SetIVModel(const G4String& v) { fIVModel = v; }

  void SetIVQuadField(G4double v)    { fIVQuadField = v; }
  void SetIVQuadRate(G4double v)     { fIVQuadRate = v; }
  void SetIVQuadExponent(G4double v) { fIVQuadExponent = v; }

  void SetIVLinRate0(G4double v)     { fIVLinRate0 = v; }
  void SetIVLinRate1(G4double v)     { fIVLinRate1 = v; }
  void SetIVLinExponent(G4double v)  { fIVLinExponent = v; }

  void SetAlpha(G4double v)	     { fAlpha = v; }
  void SetAcousticDeform(G4double v) { fAcDeform = v; }
  void SetIVDeform(const std::vector<G4double>& vlist) { fIVDeform = vlist; }
  void SetIVEnergy(const std::vector<G4double>& vlist) { fIVEnergy = vlist; }

  const G4String& GetIVModel() const { return fIVModel; }

  G4double GetIVQuadField() const    { return fIVQuadField; }
  G4double GetIVQuadRate() const     { return fIVQuadRate; }
  G4double GetIVQuadExponent() const { return fIVQuadExponent; }

  G4double GetIVLinRate0() const     { return fIVLinRate0; }
  G4double GetIVLinRate1() const     { return fIVLinRate1; }
  G4double GetIVLinExponent() const  { return fIVLinExponent; }

  G4double GetAlpha() const	     { return fAlpha; }
  G4double GetAcousticDeform() const { return fAcDeform; }
  G4int    GetNIVDeform() const { return (G4int)fIVDeform.size(); }
  const std::vector<G4double>& GetIVDeform() const { return fIVDeform; }
  const std::vector<G4double>& GetIVEnergy() const { return fIVEnergy; }
  G4double GetIVDeform(G4int i) const {
    return (i>=0 && i<GetNIVDeform()) ? fIVDeform[i] : 0.;
  }
  G4double GetIVEnergy(G4int i) const {
    return (i>=0 && i<GetNIVDeform()) ? fIVEnergy[i] : 0.;
  }

private:
  void CheckBasis();	// Initialize or complete (via cross) basis vectors
  void FillElasticity();	// Unpack reduced Cij into full Cijlk
  void FillMaps();	// Populate lookup tables using kinematics calculator
  void FillMassInfo();	// Called from SetMassTensor() to compute derived forms

  // Get theta, phi bins and offsets for interpolation
  G4bool FindLookupBins(const G4ThreeVector& k, G4int& iTheta, G4int& iPhi,
			G4double& dTheta, G4double& dPhi) const;

  // Use lookup table to get group velocity for phonons
  G4ThreeVector LookupKtoVg(G4int mode, const G4ThreeVector& k) const;

  // Use direct calculation to get group velocity for phonons
  G4ThreeVector ComputeKtoVg(G4int mode, const G4ThreeVector& k) const;

private:
  // Create a thread-local buffer to use with MapAtoB() functions
  inline G4ThreeVector& tempvec() const {
    static G4ThreadLocal G4ThreeVector* v=0;
    if (!v) v = new G4ThreeVector;
    return *v;
  }

private:
  mutable G4int verboseLevel;		    // Enable diagnostic output
  G4String fName;			    // Name of lattice for messages
  G4CMPCrystalGroup fCrystal;		    // Symmetry group, axis unit vectors
  G4ThreeVector fBasis[3];		    // Basis vectors for Miller indices
  G4double fDensity;			    // Material density (natural units)
  G4double fNImpurity;			    // Neutral impurity number density
  G4double fPermittivity;		    // Material epsilon/epsilon0 
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
  G4double fTTFrac;  // Fraction of anharmonic decays L -> TT
  G4double fBeta, fGamma, fLambda, fMu; // dynamical constants for material
  G4double fDebye;   // Debye energy, for partitioning primary phonons

  G4double fVSound;	// Speed of sound (longitudinal phonon)
  G4double fVTrans;	// Speed of sound (transverse phonon)
  G4double fL0_e;	// Scattering length for electrons
  G4double fL0_h;	// Scattering length for holes

  const G4double mElectron;	 // Free electron mass (without G4's c^2)
  G4double fHoleMass;		 // Effective mass of +ve carrier
  G4double fElectronMass;	 // Effective mass (scalar) of -ve carrier
  G4double fElectronMDOS;	 // Density of states weighed -ve carrier mass

  G4double fBandGap;	 // Minimum band gap energy
  G4double fPairEnergy;   // electron-hole pair production average energy
  G4double fFanoFactor;   // Fano factor (duh)
  G4RotationMatrix fMassTensor;	 // Full electron mass tensor
  G4RotationMatrix fMassInverse; // Inverse electron mass tensor (convenience)
  G4RotationMatrix fMassRatioSqrt;       // SQRT of tensor/scalar ratio
  G4RotationMatrix fMInvRatioSqrt;       // SQRT of scalar/tensor ratio
  std::vector<G4RotationMatrix> fValley; // Electron transport directions
  std::vector<G4RotationMatrix> fValleyInv;
  std::vector<G4ThreeVector> fValleyAxis;

  G4double fAlpha;			// Non-parabolicity of -ve potential
  G4double fAcDeform;		 	// Deformation potential for acoustic IV
  std::vector<G4double> fIVDeform;	// D0, D1 potentials for optical IV
  std::vector<G4double> fIVEnergy;	// D0, D1 thresholds for optical IV

  G4double fIVQuadField;	 // Edelweiss field scale for IV scattering
  G4double fIVQuadRate;		 // Edelweiss rate factor for IV scattering
  G4double fIVQuadExponent;	 // Edelweiss power law for E-field in IV
  G4double fIVLinExponent;	 // Power law for linear scaled IV scattering
  G4double fIVLinRate0;		 // Constant rate for linear scaled IV scat.
  G4double fIVLinRate1;		 // Linear rate for linear scaled IV scat.

  G4String fIVModel;		 // Name of IV rate function to be used
};

// Write lattice structure to output stream

inline std::ostream& 
operator<<(std::ostream& os, const G4LatticeLogical& lattice) {
  lattice.Dump(os);
  return os;
}

#endif
