/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/include/G4LatticePhysical.hh
/// \brief Definition of the G4LatticePhysical class
//
// $Id$
//
// 20131114  Add verbosity for diagnostic output
// 20131116  Replace G4Transform3D with G4RotationMatrix
// 20140312  Add pass-through calls for charge-carrier support
// 20140319  Add output functions for diagnostics
// 20140321  Move placement transformations to G4CMPProcessUtils, put
//		lattice orientation into ctor arguments
// 20140324  Add intervalley scattering parameters
// 20140401  Add valley momentum calculations
// 20140425  Add "effective mass" calculation for electrons
// 20150601  Add mapping from electron velocity back to momentum
// 20160608  Drop (theta,phi) lattice orientation function.
// 20170523  Add interface for axis vector of valleys
// 20170525  Drop empty destructor to allow default "rule of five" semantics
// 20170810  Add call-throughs for new IV scattering parameters
// 20170821  Add transverse sound speed, L->TT fraction
// 20170919  Add access to full lists of IV scattering matrix terms
// 20170928  Replace "pol" with "mode" for phonon states
// 20180831  Add accessors for linear IV rate, change Edelweiss names
// 20181001  M. Kelsey -- Clarify IV rate parameters systematically
// 20190704  M. Kelsey -- Add IV rate function selector for material
// 20190801  M. Kelsey -- Use G4ThreeVector buffer instead of pass-by-value;
//		add pass through call to GetValleyInv().
// 20200520  For MT thread safety, wrap G4ThreeVector buffer in function to
//		return thread-local instance.
// 20200608  Fix -Wshadow warnings from tempvec
// 20210919  M. Kelsey -- Allow SetVerboseLevel() from const instances.
// 20220921  G4CMP-319 -- Add utilities for thermal (Maxwellian) distributions
//		Also, add long missing accessors for Miller orientation

#ifndef G4LatticePhysical_h
#define G4LatticePhysical_h 1

#include "G4LatticeLogical.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>

#define G4CMP_HAS_TEMPERATURE	/* G4CMP-319 -- New feature for user code */


class G4LatticePhysical {
public:
  G4LatticePhysical();		// User *MUST* set configuration manually

  // Miller orientation aligns lattice normal (hkl) with geometry +Z
  G4LatticePhysical(const G4LatticeLogical* Lat,
		    G4int h=0, G4int k=0, G4int l=0, G4double rot=0.);

  void SetVerboseLevel(G4int vb) const {
    verboseLevel = vb;
    if (fLattice) fLattice->SetVerboseLevel(vb);
  }

  // Specific material lattice for this physical instance
  void SetLatticeLogical(const G4LatticeLogical* Lat) { fLattice = Lat; }

  // Set physical lattice orientation, relative to G4VSolid coordinates
  // Miller orientation aligns lattice normal (hkl) with geometry +Z
  void SetMillerOrientation(G4int h, G4int k, G4int l, G4double rot=0.);

  // Set temperature of volume/lattice for use with thermalization processes
  void SetTemperature(G4double temp) { fTemperature = temp; }

  // Rotate input vector between lattice and solid orientations
  // Returns new vector value for convenience
  const G4ThreeVector& RotateToLattice(G4ThreeVector& dir) const;
  const G4ThreeVector& RotateToSolid(G4ThreeVector& dir) const;

  // Convert input wave vector and polarization to group velocity
  // NOTE:  Input vector must be in local (G4VSolid) coordinate system
  G4double      MapKtoV(G4int mode, const G4ThreeVector& k) const;
  G4ThreeVector MapKtoVDir(G4int mode, const G4ThreeVector& k) const;

  // Convert between electron momentum and valley velocity or HV wavevector
  // NOTE:  p or v_el vector must be in local (G4VSolid) coordinate system
  // NOTE:  K_HV vector must be in valley internal coordinate system
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

public:  
  const G4LatticeLogical* GetLattice() const { return fLattice; }

  // Return Miller orientation
  G4double GetRotation() const { return fRot; }
  G4ThreeVector GetMillerIndices() const {
    return G4ThreeVector(hMiller, kMiller, lMiller);
  }

  // Return temperature assigned to lattice/volume, or global setting
  G4double GetTemperature() const;

  // Call through to get material properties
  G4double GetDensity() const { return fLattice->GetDensity(); }
  G4double GetImpurities() const { return fLattice->GetImpurities(); }
  G4double GetPermittivity() const { return fLattice->GetPermittivity(); }

  // Call through to get crystal basis vectors
  const G4ThreeVector& GetBasis(G4int i) const { return fLattice->GetBasis(i); }

  // Phonon propagation parameters
  G4double GetScatteringConstant() const { return fLattice->GetScatteringConstant(); }
  G4double GetAnhDecConstant() const { return fLattice->GetAnhDecConstant(); }
  G4double GetAnhTTFrac() const      { return fLattice->GetAnhTTFrac(); }
  G4double GetLDOS() const           { return fLattice->GetLDOS(); }
  G4double GetSTDOS() const          { return fLattice->GetSTDOS(); }
  G4double GetFTDOS() const          { return fLattice->GetFTDOS(); }
  G4double GetBeta() const           { return fLattice->GetBeta(); }
  G4double GetGamma() const          { return fLattice->GetGamma(); }
  G4double GetLambda() const         { return fLattice->GetLambda(); }
  G4double GetMu() const             { return fLattice->GetMu(); }
  G4double GetDebyeEnergy() const    { return fLattice->GetDebyeEnergy(); }

  // Charge carrier propagation parameters
  G4double GetBandGapEnergy() const   { return fLattice->GetBandGapEnergy(); }
  G4double GetPairProductionEnergy() const { return fLattice->GetPairProductionEnergy(); }
  G4double GetFanoFactor() const      { return fLattice->GetFanoFactor(); }
  G4double GetSoundSpeed() const      { return fLattice->GetSoundSpeed(); }
  G4double GetTransverseSoundSpeed() const { return fLattice->GetTransverseSoundSpeed(); }
  G4double GetElectronScatter() const { return fLattice->GetElectronScatter(); }
  G4double GetHoleScatter() const     { return fLattice->GetHoleScatter(); }

  // Charge carriers have effective mass
  G4double GetHoleMass() const { return fLattice->GetHoleMass(); }
  G4double GetElectronMass() const { return fLattice->GetElectronMass(); }
  G4double GetElectronDOSMass() const { return fLattice->GetElectronDOSMass(); }
  G4double GetElectronEffectiveMass(G4int iv, const G4ThreeVector& p) const {
    return fLattice->GetElectronEffectiveMass(iv, p);
  }

  const G4RotationMatrix& GetMassTensor() const { return fLattice->GetMassTensor(); }
  const G4RotationMatrix& GetMInvTensor() const { return fLattice->GetMInvTensor(); }
  const G4RotationMatrix& GetSqrtTensor() const { return fLattice->GetSqrtTensor(); }
  const G4RotationMatrix& GetSqrtInvTensor() const { return fLattice->GetSqrtInvTensor(); }

  // Electrons are biased to move along energy minima in momentum space
  size_t NumberOfValleys() const { return fLattice->NumberOfValleys(); }

  // FIXME:  Should valley matrix be rotated from internal to local coordinates?
  const G4RotationMatrix& GetValley(G4int iv) const    { return fLattice->GetValley(iv); }
  const G4RotationMatrix& GetValleyInv(G4int iv) const { return fLattice->GetValleyInv(iv); }
  const G4ThreeVector& GetValleyAxis(G4int iv) const   { return fLattice->GetValleyAxis(iv); }

  // Parameters for electron intervalley scattering (Edelweiss, linear, matrix)
  const G4String& GetIVModel() const { return fLattice->GetIVModel(); }

  G4double GetIVQuadField() const    { return fLattice->GetIVQuadField(); }
  G4double GetIVQuadRate() const     { return fLattice->GetIVQuadRate(); }
  G4double GetIVQuadExponent() const { return fLattice->GetIVQuadExponent(); }

  G4double GetIVLinRate0() const    { return fLattice->GetIVLinRate0(); }
  G4double GetIVLinRate1() const    { return fLattice->GetIVLinRate1(); }
  G4double GetIVLinExponent() const { return fLattice->GetIVLinExponent(); }

  G4double GetAlpha() const          { return fLattice->GetAlpha(); }
  G4double GetAcousticDeform() const { return fLattice->GetAcousticDeform(); }

  // Optical intervalley scattering may use D0 or D1 deformation potentials
  G4int    GetNIVDeform() const { return fLattice->GetNIVDeform(); }
  G4double GetIVDeform(G4int i) const { return fLattice->GetIVDeform(i); }
  G4double GetIVEnergy(G4int i) const { return fLattice->GetIVEnergy(i); }
  const std::vector<G4double>& GetIVDeform() const { return fLattice->GetIVDeform(); }
  const std::vector<G4double>& GetIVEnergy() const { return fLattice->GetIVEnergy(); }

  // Dump logical lattice, with additional info about physical
  void Dump(std::ostream& os) const;

private:
  // Create a thread-local buffer to use with MapAtoB() functions
  inline G4ThreeVector& tempvec() const {
    static G4ThreadLocal G4ThreeVector* v=0;
    if (!v) v = new G4ThreeVector;
    return *v;
  }

private:
  mutable G4int verboseLevel;		// Enable diagnostic output
  const G4LatticeLogical* fLattice;	// Underlying lattice parameters
  G4RotationMatrix fOrient;		// Rotate geometry into lattice frame
  G4RotationMatrix fInverse;
  G4int hMiller, kMiller, lMiller;	// Save Miller indices for dumps
  G4double fRot;
  G4double fTemperature;		// Temperature assigned to volume
};

// Write lattice structure to output stream

inline std::ostream& 
operator<<(std::ostream& os, const G4LatticePhysical& lattice) {
  lattice.Dump(os);
  return os;
}

#endif
