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

#ifndef G4LatticePhysical_h
#define G4LatticePhysical_h 1

#include "G4LatticeLogical.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>


class G4LatticePhysical {
public:
  G4LatticePhysical();		// User *MUST* set configuration manually

  G4LatticePhysical(const G4LatticeLogical* Lat,	// Lattice orientation
		    G4double theta=0., G4double phi=0.);

  G4LatticePhysical(const G4LatticeLogical* Lat,	// Miller indices
		    G4int h, G4int k, G4int l);

  virtual ~G4LatticePhysical();

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  // Specific material lattice for this physical instance
  void SetLatticeLogical(const G4LatticeLogical* Lat) { fLattice = Lat; }

  // Set physical lattice orientation, relative to G4VSolid coordinates
  void SetLatticeOrientation(G4double theta, G4double phi);
  void SetMillerOrientation(G4int h, G4int k, G4int l);

  // Rotate input vector between lattice and solid orientations
  // Returns new vector value for convenience
  const G4ThreeVector& RotateToLattice(G4ThreeVector& dir) const;
  const G4ThreeVector& RotateToSolid(G4ThreeVector& dir) const;

  // Convert input wave vector and polarization to group velocity
  // NOTE:  Input vector must be in local (G4VSolid) coordinate system
  // NOTE:  Pass vector by value to allow in-situ rotations
  G4double      MapKtoV(G4int pol, G4ThreeVector k) const;
  G4ThreeVector MapKtoVDir(G4int pol, G4ThreeVector k) const;

  // Convert between electron momentum and valley velocity or HV wavevector
  // NOTE:  Input vector must be in local (G4VSolid) coordinate system
  // NOTE:  Pass vector by value to allow in-situ rotations
  G4ThreeVector MapPtoV_el(G4int ivalley, G4ThreeVector p_e) const;
  G4ThreeVector MapV_elToP(G4int ivalley, G4ThreeVector v_el) const;
  G4ThreeVector MapV_elToK_HV(G4int ivalley, G4ThreeVector v_el) const;
  G4ThreeVector MapPtoK_valley(G4int ivalley, G4ThreeVector p_e) const;
  G4ThreeVector MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const;
  G4ThreeVector MapK_HVtoP(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_HVtoK_valley(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_HVtoK(G4int ivalley, G4ThreeVector k_HV) const;
  G4ThreeVector MapK_valleyToP(G4int ivalley, G4ThreeVector k) const;

  // Apply energy relationships for electron transport
  G4double MapPtoEkin(G4int ivalley, G4ThreeVector p_e) const;
  G4double MapV_elToEkin(G4int ivalley, G4ThreeVector v_e) const;

public:  
  const G4LatticeLogical* GetLattice() const { return fLattice; }

  // Call through to get crystal basis vectors
  const G4ThreeVector& GetBasis(G4int i) const { return fLattice->GetBasis(i); }

  // Phonon propagation parameters
  G4double GetScatteringConstant() const { return fLattice->GetScatteringConstant(); }
  G4double GetAnhDecConstant() const { return fLattice->GetAnhDecConstant(); }
  G4double GetLDOS() const           { return fLattice->GetLDOS(); }
  G4double GetSTDOS() const          { return fLattice->GetSTDOS(); }
  G4double GetFTDOS() const          { return fLattice->GetFTDOS(); }
  G4double GetBeta() const           { return fLattice->GetBeta(); }
  G4double GetGamma() const          { return fLattice->GetGamma(); }
  G4double GetLambda() const         { return fLattice->GetLambda(); }
  G4double GetMu() const             { return fLattice->GetMu(); }

  // Charge carrier propagation parameters
  G4double GetSoundSpeed() const      { return fLattice->GetSoundSpeed(); }
  G4double GetElectronScatter() const { return fLattice->GetElectronScatter(); }
  G4double GetHoleScatter() const     { return fLattice->GetHoleScatter(); }

  // Charge carriers have effective mass
  G4double GetHoleMass() const { return fLattice->GetHoleMass(); }
  G4double GetElectronMass() const { return fLattice->GetElectronMass(); }
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
  const G4RotationMatrix& GetValley(G4int iv) const { return fLattice->GetValley(iv); }

  // Parameters for electron intervalley scattering
  G4double GetIVField() const    { return fLattice->GetIVField(); }
  G4double GetIVRate() const     { return fLattice->GetIVRate(); }
  G4double GetIVExponent() const { return fLattice->GetIVExponent(); }

  // Dump logical lattice, with additional info about physical
  void Dump(std::ostream& os) const;

private:
  G4int verboseLevel;			// Enable diagnostic output

  G4double fTheta, fPhi;		// Lattice orientation within object
  const G4LatticeLogical* fLattice;	// Underlying lattice parameters
};

// Write lattice structure to output stream

inline std::ostream& 
operator<<(std::ostream& os, const G4LatticePhysical& lattice) {
  lattice.Dump(os);
  return os;
}

#endif
