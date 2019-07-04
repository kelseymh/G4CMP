/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170822 M. Kelsey -- Rename EnergyLimiter to TrackLimiter

#ifndef G4CMPProcessSubType_hh
#define G4CMPProcessSubType_hh 1


// NOTE:  SubType codes have to be globally unique; bad design!

enum G4CMPProcessSubType {
  fPhononScattering = 301,
  fPhononReflection,
  fPhononDownconversion,
  fInterValleyScattering,
  fLukeScattering,
  fChargeBoundary,
  fTimeStepper,
  fSecondaryProduction,
  fTrackLimiter,
  fChargeRecombine
};

#endif	/* G4CMPProcessSubType_hh */
