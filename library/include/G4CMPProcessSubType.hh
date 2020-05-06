/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170822 M. Kelsey -- Rename EnergyLimiter to TrackLimiter
// 20200331 C. Stanford G4CMP-195:  Add Trapping and Impact subtypes
// 20200504 M. Kelsy -- Remove impact subtype here; set values explicitly

#ifndef G4CMPProcessSubType_hh
#define G4CMPProcessSubType_hh 1


// NOTE:  SubType codes have to be globally unique; bad design!

enum G4CMPProcessSubType {
  fPhononScattering = 301,
  fPhononReflection = 302,
  fPhononDownconversion = 303,
  fInterValleyScattering = 304,
  fLukeScattering = 305,
  fChargeBoundary = 306,
  fTimeStepper = 307,
  fSecondaryProduction = 308,
  fTrackLimiter = 309,
  fChargeRecombine = 310,
  fChargeTrapping = 313
};

#endif	/* G4CMPProcessSubType_hh */
