/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPIVRateQuadratic.cc
/// \brief Compute electron intervalley scattering rate using linear
///	   power-law parametrization vs. electric field.
//
// $Id$
//
// 20181001  Use systematic names for IV rate parameters
// 20210908  Use global track position to query field; configure field.
// 20211003  Use encapsulated G4CMPFieldUtils to get field.
// 20230829  Rotated E-field to local frame first and changed Mass Multiplication in HV transformation

#include "G4CMPIVRateLinear.hh"
#include "G4CMPFieldUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include <math.h>
#include <iostream>

// Scattering rate is computed from electric field

G4double G4CMPIVRateLinear::Rate(const G4Track& aTrack) const {
  G4ThreeVector fieldVector = G4CMP::GetFieldAtPosition(aTrack);

  if (verboseLevel > 1) {
    G4cout << "IV local position " << GetLocalPosition(aTrack) << G4endl
	   << " field " << fieldVector/volt*cm << " V/cm" << G4endl
	   << " magnitude " << fieldVector.mag()/volt*cm << " V/cm toward "
	   << fieldVector.cosTheta() << " z" << G4endl;
  }

  // Find E-field in HV space: in lattice frame, rotate into valley,
  // then apply HV tansform.
  // NOTE:  Separate steps to avoid matrix-matrix multiplications
  fieldVector = theLattice->EllipsoidalToSphericalTranformation(GetValleyIndex(aTrack), fieldVector);
  // fieldVector *= sqrt(theLattice->GetElectronMass()/(electron_mass_c2/c_squared));
  fieldVector /= volt/cm;			// Strip units for MFP below
  if (verboseLevel > 1) {
    G4cout << " in HV space " << fieldVector << " ("
	   << fieldVector.mag() << ") V/cm" << G4endl;
  }

  // Compute mean free path -- NOTE FIELD UNITS ARE V/cm HERE
  G4double rate = theLattice->GetIVLinRate0() +
    theLattice->GetIVLinRate1() * pow(fieldVector.mag(),
				      theLattice->GetIVLinExponent());

  if (verboseLevel > 1) G4cout << "IV rate = " << rate/hertz << " Hz" << G4endl;
  return rate;
}
