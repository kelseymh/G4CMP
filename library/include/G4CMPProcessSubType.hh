/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20170822 M. Kelsey -- Rename EnergyLimiter to TrackLimiter
// 20200501 G4CMP-196: Need separate processes for A- and D- charge traps

#ifndef G4CMPProcessSubType_hh
#define G4CMPProcessSubType_hh 1


// NOTE 1:  SubType codes have to be globally unique; bad design!
// NOTE 2:  Enumerator valus given explicitly below to help w/OrdParamTable

enum G4CMPProcessSubType {
  fG4CMPProcess = 300,		// Use this like "unknown", not specific
  fPhononScattering = 301,
  fPhononReflection = 302,
  fPhononDownconversion = 303,
  fInterValleyScattering =304,
  fLukeScattering = 305,
  fChargeBoundary = 306,
  fTimeStepper = 307,
  fSecondaryProduction = 308,
  fTrackLimiter = 309,
  fChargeRecombine = 310,
  fDTrapIonization = 311,
  fATrapIonization = 312,
};

#endif	/* G4CMPProcessSubType_hh */
