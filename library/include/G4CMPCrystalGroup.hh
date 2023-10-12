#ifndef G4CMPCrystalGroup_hh
#define G4CMPCrystalGroup_hh 1
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
/// \file library/include/G4CMPCrystalGroup.hh
/// \brief Enumerator and support functions for lattice symmetry groups
//
// 20160729  M. Kelsey -- Add accessors for unit cell angles
// 20170525  M. Kelsey -- Add default "rule of five" copy/move operators
// 20170728  Change function args "alpha, beta, gamma" to "al, bt, gm" (-Wshadow)

#include "globals.hh"
#include "G4ThreeVector.hh"


class G4CMPCrystalGroup {
public:
  enum Bravais { amorphous, cubic, tetragonal, orthorhombic, hexagonal,
		 rhombohedral, monoclinic, triclinic, UNKNOWN=-1 };
  static const char* Name(Bravais grp);
  static Bravais Group(const G4String& name);	// May return -1 if invalid

public:				// For convenient access to data members 
  Bravais group;
  G4ThreeVector axis[3];	// Basis unit vectors in direct orientation

public:
  G4CMPCrystalGroup() : group(UNKNOWN) {;}	// Default ctor, must use Set()
  G4CMPCrystalGroup(Bravais grp) { Set(grp); }
  G4CMPCrystalGroup(Bravais grp, G4double angle) { Set(grp, angle); }
  G4CMPCrystalGroup(G4double al, G4double bt, G4double gm) {
    Set(triclinic, al, bt, gm);
  }

  virtual ~G4CMPCrystalGroup() {;}

  G4CMPCrystalGroup(const G4CMPCrystalGroup&) = default;
  G4CMPCrystalGroup(G4CMPCrystalGroup&&) = default;
  G4CMPCrystalGroup& operator=(const G4CMPCrystalGroup&) = default;
  G4CMPCrystalGroup& operator=(G4CMPCrystalGroup&&) = default;

  const char* Name() const { return Name(group); }

  G4double alpha() const { return fabs(axis[0].angle(axis[2])); }
  G4double beta() const  { return fabs(axis[1].angle(axis[2])); }
  G4double gamma() const { return fabs(axis[0].angle(axis[1])); }

  // Some parameters may be omitted depending on symmetry
  void Set(Bravais grp, G4double a=0., G4double b=0., G4double g=0.);

  // Copy appropriate elements of Cij matrix based on crystal symmetry
  // NOTE:  Non-const array passed in for modification
  G4bool FillElReduced(G4double Cij[6][6]) const;

private:
  void SetCartesian();
  void SetHexagonal();
  void SetRhombohedral(G4double angle);
  void SetMonoclinic(G4double angle);
  void SetTriclinic(G4double al, G4double bt, G4double gm);

  // Separate functions for completing Cij to avoid long spaghetti code
  G4bool FillCubic(G4double Cij[6][6]) const;
  G4bool FillTetragonal(G4double Cij[6][6]) const;
  G4bool FillOrthorhombic(G4double Cij[6][6]) const;
  G4bool FillRhombohedral(G4double Cij[6][6]) const;
  G4bool FillMonoclinic(G4double Cij[6][6]) const;
  G4bool FillTriclinic(G4double Cij[6][6]) const;
  void ReflectElReduced(G4double Cij[6][6]) const;	// Applies symmetry
};

#endif	/* G4CMPCrystalGroup_hh */
