/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPTrackLimiter.hh
/// \brief Definition of the G4CMPTrackLimiter process, to kill tracks
///        falling below a minimum energy (set in G4CMPConfigManager),
///	   or which travel outside their original volume.
//
// $Id$
//
// 20170602  M. Kelsey -- Inherit from new G4CMPVProcess
// 20170822  M. Kelsey -- Add checking on current vs. original volume

#ifndef G4CMPTrackLimiter_hh
#define G4CMPTrackLimiter_hh 1

#include "G4CMPVProcess.hh"

class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VParticleChange;


class G4CMPTrackLimiter : public G4CMPVProcess {
public:
  G4CMPTrackLimiter(const G4String& name="TrackLimiter")
    : G4CMPVProcess(name, fTrackLimiter) {;}
  virtual ~G4CMPTrackLimiter() {;}

  virtual G4bool IsApplicable(const G4ParticleDefinition& pd);

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
				       G4ForceCondition*);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  G4bool BelowEnergyCut(const G4Track& track) const;

  G4bool EscapedFromVolume(const G4Step& step) const;

  virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*);

private:
  G4CMPTrackLimiter(const G4CMPTrackLimiter&);	// Copying is forbidden
  G4CMPTrackLimiter& operator=(const G4CMPTrackLimiter&);
};

#endif	/* G4CMPTrackLimiter_hh */
