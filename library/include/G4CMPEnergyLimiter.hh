/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPEnergyLimiter.hh
/// \brief Definition of the G4CMPEnergyLimiter process, to kill tracks
///        falling below a minimum energy (set in G4CMPConfigManager).
//
// $Id$
//

#ifndef G4CMPEnergyLimiter_hh
#define G4CMPEnergyLimiter_hh 1

#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPProcessSubType.hh"


class G4CMPEnergyLimiter : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4CMPEnergyLimiter(const G4String& name="EnergyLimiter")
    : G4VDiscreteProcess(name, fPhonon) { SetProcessSubType(fEnergyLimiter); }
  virtual ~G4CMPEnergyLimiter() {;}

  virtual G4bool IsApplicable(const G4ParticleDefinition& pd);

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
				       G4ForceCondition*);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  G4bool BelowEnergyCut(const G4Track& track) const;

  virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*) {
    return DBL_MAX;
  }

private:
  G4CMPEnergyLimiter(const G4CMPEnergyLimiter&);	// Copying is forbidden
  G4CMPEnergyLimiter& operator=(const G4CMPEnergyLimiter&);
};

#endif	/* G4CMPEnergyLimiter_hh */
