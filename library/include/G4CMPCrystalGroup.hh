#ifndef G4CMPCrystalGroup_hh
#define G4CMPCrystalGroup_hh 1
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononPolarization.hh
/// \brief Enumerator and support functions for lattice symmetry groups
//
// $Id$

#include "globals.hh"
#include "G4ThreeVector.hh"


class G4CMPCrystalGroup {
public:
  enum Bravais { amorphous, cubic, tetragonal, orthorhombic, hexagonal,
		 rhombohedral, monoclinic, triclinic };
  static const char* Name(Bravais grp);

public:				// For convenient access to data members 
  Bravais group;
  G4ThreeVector axis[3];	// Basis unit vectors in direct orientation

public:
  G4CMPCrystalGroup(Bravais grp);			// Orthogonal axes only
  G4CMPCrystalGroup(Bravais grp, G4double angle);	// One non-ortho angle
  G4CMPCrystalGroup(G4double alpha, G4double beta,	// Triclinic only
		    G4double gamma);

  virtual ~G4CMPCrystalGroup() {;}

  const char* Name() const { return Name(group); }

  // Copy appropriate elements of Cij matrix based on crystal symmetry
  // NOTE:  Non-const array passed in for modification
  G4bool FillElReduced(G4double Cij[6][6]) const;

private:
  void SetCartesian();
  void SetHexagonal();
  void SetRhombohedral(G4double angle);
  void SetMonoclinic(G4double angle);
  void SetTriclinic(G4double alpha, G4double beta, G4double gamma);

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
