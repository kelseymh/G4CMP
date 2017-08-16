/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4PhononPolarization.hh
/// \brief Enumerator and support functions for lattice symmetry groups
//
// $Id$
//
// 20170728  Change function args "alpha, beta, gamma" to "al, bt, gm" (-Wshadow)

#include "G4CMPCrystalGroup.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


// Some input angles may be ignored, depending on crystal symmetry

void G4CMPCrystalGroup::Set(G4CMPCrystalGroup::Bravais grp,
			    G4double al, G4double bt, G4double gm) {
  group = grp;
  switch (group) {
  case amorphous:
  case cubic:
  case tetragonal:
  case orthorhombic: SetCartesian(); break;
  case hexagonal:    SetHexagonal(); break;
  case monoclinic:   SetMonoclinic(al); break;
  case rhombohedral: SetRhombohedral(al); break;
  case triclinic:    SetTriclinic(al, bt, gm); break;
  default:
    G4ExceptionDescription msg;
    msg << "Unrecognized crystal group " << group;
    G4Exception("G4CMPCrystalGroup", "Lattice100", FatalException, msg);
  }
}


// Fill Bravais axis unit vectors according to crystal system

void G4CMPCrystalGroup::SetCartesian() {
  axis[0] = CLHEP::HepXHat;
  axis[1] = CLHEP::HepYHat;
  axis[2] = CLHEP::HepZHat;
}

void G4CMPCrystalGroup::SetHexagonal() {
  SetCartesian();
  axis[1].rotateZ(30.*deg);		// X-Y opening angle is 120 degrees
}

void G4CMPCrystalGroup::SetRhombohedral(G4double angle) {
  SetTriclinic(angle, angle, angle);
}
 
void G4CMPCrystalGroup::SetMonoclinic(G4double angle) {
  SetCartesian();
  axis[2].rotateX(angle-halfpi);	// Z-Y opening angle
}
 
void G4CMPCrystalGroup::SetTriclinic(G4double al, G4double bt,
				     G4double gm) {
  SetCartesian();
  axis[1].rotateZ(gm-halfpi);	// X-Y opening angle

  // Z' axis computed by hand to get both opening angles right
  // X'.Z' = cos(alpha), Y'.Z' = cos(beta), solve for Z' components
  G4double ca=cos(al), cb=cos(bt), cg=cos(gm), sg=sin(gm);
  G4double x3=ca, y3=(cb-ca*cg)/sg, z3=sqrt(1.-x3*x3-y3*y3);
  axis[2].set(x3, y3, z3);
}


// Copy appropriate elements of Cij matrix based on crystal symmetry
// Returns false if some expected element is zero.

G4bool G4CMPCrystalGroup::FillElReduced(G4double Cij[6][6]) const {
  switch (group) {
  case amorphous:    Cij[3][3] = 0.5*(Cij[0][0]-Cij[0][1]); // Cubic, C44 set
  case cubic:        return FillCubic(Cij); break;
  case hexagonal:    Cij[0][5] = 0.;			// Tetragonal, C16=0
                     Cij[4][5] = 0.5*(Cij[0][0] - Cij[0][1]);
  case tetragonal:   return FillTetragonal(Cij); break;
  case orthorhombic: return FillOrthorhombic(Cij); break;
  case rhombohedral: return FillRhombohedral(Cij); break;
  case monoclinic:   return FillMonoclinic(Cij); break;
  case triclinic:    return FillTriclinic(Cij); break;
  default: break;
  }

  return false;
}

G4bool G4CMPCrystalGroup::FillCubic(G4double Cij[6][6]) const {
  G4double C11=Cij[0][0], C12=Cij[0][1], C44=Cij[3][3];

  for (size_t i=0; i<6; i++) {
    for (size_t j=i; j<6; j++) {
      if (i<3 && j<3) Cij[i][j] = (i==j) ? C11 : C12;
      else if (i==j && i>=3) Cij[i][i] = C44;	    
      else Cij[i][j] = 0.;
    }
  }

  ReflectElReduced(Cij);

  return (C11!=0. && C12!=0. && C44!=0.);
}

G4bool G4CMPCrystalGroup::FillTetragonal(G4double Cij[6][6]) const {
  G4double C11=Cij[0][0], C12=Cij[0][1], C13=Cij[0][2], C16=Cij[0][5];
  G4double C33=Cij[2][2], C44=Cij[3][3], C66=Cij[5][5];

  Cij[1][1] = C11;	// Copy small number of individual elements
  Cij[1][2] = C13;
  Cij[1][5] = -C16;
  Cij[4][4] = C44;

  ReflectElReduced(Cij);

  // NOTE:  Do not test for C16 != 0., to allow calling from Hexagonal
  return (C11!=0. && C12!=0. && C13!=0. && C33!=0. && C44!=0. && C66!=0.);
}

G4bool G4CMPCrystalGroup::FillOrthorhombic(G4double Cij[6][6]) const {
  // No degenerate elements; just check for all non-zero
  ReflectElReduced(Cij);

  G4bool good = true;
  for (size_t i=0; i<6; i++) {
    for (size_t j=i+1; j<3; j++) good &= (Cij[i][j] != 0);
    good &= (Cij[i][i] != 0);
  }

  return good;
}

G4bool G4CMPCrystalGroup::FillRhombohedral(G4double Cij[6][6]) const {
  G4double C11=Cij[0][0], C12=Cij[0][1], C13=Cij[0][2], C14=Cij[0][3];
  G4double C15=Cij[0][4], C33=Cij[2][2], C44=Cij[3][3], C66=0.5*(C11-C12);

  Cij[1][1] = C11;	// Copy small number of individual elements
  Cij[1][2] = C13;
  Cij[1][3] = -C14;
  Cij[1][4] = -C15;
  Cij[3][5] = -C15;
  Cij[4][4] = C44;
  Cij[4][5] = C14;

  // NOTE:  C15 may be zero (c.f. rhombohedral(I) vs. (II))
  return (C11!=0 && C12!=0 && C13!=0 && C14!=0. &&
	  C33!=0. && C44!=0. && C66!=0.);
}

G4bool G4CMPCrystalGroup::FillMonoclinic(G4double Cij[6][6]) const {
  // The monoclinic matrix has 13 independent elements with no degeneracies
  // Sanity condition is same as orthorhombic, plus C45, C(1,2,3)6

  return (FillOrthorhombic(Cij) && Cij[0][5]!=0. && Cij[1][5]!=0. &&
	  Cij[2][5] != 0. && Cij[3][4]!=0.);
}

G4bool G4CMPCrystalGroup::FillTriclinic(G4double Cij[6][6]) const {
  // The triclinic matrix has the entire upper half filled (21 elements)

  ReflectElReduced(Cij);

  G4bool good = true;
  for (size_t i=0; i<6; i++) {
    for (size_t j=i; j<6; j++) good &= (Cij[i][j] != 0);
  }

  return good;
}


// Apply matrix symmetry about diagonal: Cji = Cij
void G4CMPCrystalGroup::ReflectElReduced(G4double Cij[6][6]) const {
  for (size_t i=0; i<6; i++) {
    for (size_t j=i+1; j<6; j++) {
      Cij[j][i] = Cij[i][j];
    }
  }
}


// Convert between enumerator and useful strings

const char* G4CMPCrystalGroup::Name(Bravais grp) {
  switch (grp) {
  case amorphous:    return "amorphous";
  case cubic:        return "cubic";
  case tetragonal:   return "tetragonal";
  case orthorhombic: return "orthorhombic";
  case hexagonal:    return "hexagonal";
  case rhombohedral: return "rhombohedral";
  case monoclinic:   return "monoclinic";
  case triclinic:    return "triclinic";
  default: break;
  }

  return 0;	// Force crash: should never get here if passed enum
}

G4CMPCrystalGroup::Bravais G4CMPCrystalGroup::Group(const G4String& name) {
  if (name.index("amo")==0) return amorphous;
  if (name.index("cub")==0) return cubic;
  if (name.index("tet")==0) return tetragonal;
  if (name.index("ort")==0) return orthorhombic;
  if (name.index("hex")==0) return hexagonal;
  if (name.index("mon")==0) return monoclinic;
  if (name.index("tri")==0) return triclinic;

  return UNKNOWN;	// Failure condition; calling code should test
}
