/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSecondaryUtils.hh
/// \brief Free standing helper functions for creating secondary tracks
///
//
// $Id$
//
// 20161115 Initial commit - R. Agnese
// 20170620 M. Kelsey -- Replace PV arg with Touchable, for transforms
// 20170815 M. Kelsey -- Move AdjustSecondaryPosition to GeometryUtils

#ifndef G4CMPSecondaryUtils_hh
#define G4CMPSecondaryUtils_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleDefinition;
class G4Track;
class G4VTouchable;


namespace G4CMP {
  G4Track* CreateSecondary(const G4Track& track, G4ParticleDefinition* pd,
			   const G4ThreeVector& waveVec,
			   G4double energy);

  G4Track* CreatePhonon(const G4VTouchable* touch, G4int polarization,
			const G4ThreeVector& waveVec, G4double energy,
			G4double time, const G4ThreeVector& pos);
  
  G4Track* CreateChargeCarrier(const G4VTouchable* touch, G4int charge,
			       G4int valley, G4double Ekin, G4double time,
			       const G4ThreeVector& pdir,
			       const G4ThreeVector& pos);

  G4Track* CreateChargeCarrier(const G4VTouchable* touch, G4int charge,
			       G4int valley, G4double time,
			       const G4ThreeVector& p,
			       const G4ThreeVector& pos);
}

#endif	/* G4CMPSecondaryUtils_hh */
