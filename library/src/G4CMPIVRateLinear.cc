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
// 20191106  Replace field management stuff with GetFieldAtPosition().

#include "G4CMPIVRateLinear.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>
#include <iostream>

// Scattering rate is computed from electric field

G4double G4CMPIVRateLinear::Rate(const G4Track& aTrack) const {
  // Get electric field associated with current volume, if any
  G4ThreeVector fieldVector = G4CMP::GetFieldAtPosition(aTrack);

  if (verboseLevel > 1) {
    G4cout << "IV global position " << aTrack.GetPosition()
	   << "\n field " << fieldVector/volt*cm << " V/cm"
	   << "\n magnitude " << fieldVector.mag()/volt*cm << " V/cm toward "
	   << fieldVector.cosTheta() << " z" << G4endl;
  }

  // Find E-field in HV space: in lattice frame, rotate into valley,
  // then apply HV tansform.
  // NOTE:  Separate steps to avoid matrix-matrix multiplications
  theLattice->RotateToLattice(fieldVector);
  fieldVector *= GetValley(aTrack);
  fieldVector *= theLattice->GetSqrtInvTensor();
  fieldVector /= volt/cm;			// Strip units for MFP below
  if (verboseLevel > 1) {
    G4cout << " in HV space " << fieldVector << " ("
	   << fieldVector.mag() << ") V/cm" << G4endl;
  }

  // Compute mean free path -- NOTE FIELD UNITS ARE V/cm HERE
  G4double rate = theLattice->GetIVLinRate0() +
    theLattice->GetIVLinRate1() * pow(fieldVector.mag(),
				      theLattice->GetIVLinExponent());

  G4double maxRate = theLattice->GetIVMaxRate();

  if (verboseLevel > 1) {
    G4cout << "IV rate = " << rate/hertz << " Hz"
	   << " vs. max " << maxRate/hertz << G4endl;
  }

  return std::min(rate, maxRate);
}
