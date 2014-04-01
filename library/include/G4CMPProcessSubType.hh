// $Id$

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
  fTimeStepper
};

#endif	/* G4CMPProcessSubType_hh */
