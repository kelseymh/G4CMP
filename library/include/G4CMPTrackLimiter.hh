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
// 20250421  M. Kelsey -- Add comparison of track position with volume.
// 20250501  G4CMP-358 -- Identify and stop charge tracks stuck in field.
// 20250506  Add local caches to compute cumulative flight distance, RMS
// 20250515  Add configuration settings for "stuck track" cuts
// 20250801  G4CMP-326:  Kill thermal phonons if finite temperature set.
// 20251025  G4CMP-520:  Remove redundant (and incorrect) InvalidPosition().

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
    : G4CMPVProcess(name, fTrackLimiter), stepWindow(10000), minPosShift(0.),
      minFlightRMS(0.), maxPathScale(20.), flightAvg(-1.), flightAvg2(-1.),
      lastFlight(-1.), lastRMS(-1.) {;}

  virtual ~G4CMPTrackLimiter() {;}

  virtual G4bool IsApplicable(const G4ParticleDefinition& pd);

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
				       G4ForceCondition*);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Initialize flight distance caches for stuck-track evaluation
  virtual void LoadDataForTrack(const G4Track* track,
				const G4bool overrideMomentumReset=false);

protected:
  G4bool BelowEnergyCut(const G4Track& track) const;
  G4bool EscapedFromVolume(const G4Step& step) const;
  G4bool ChargeStuck(const G4Track& track);	// Non-const to use caches
  G4bool PhononIsThermal(const G4Track& track) const;

  virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*);

protected:
  G4int stepWindow;		// Number of steps to use for averaging
  G4double minPosShift;		// Minimum allowed position shift by track
  G4double minFlightRMS;	// Minimum allowed RMS of flight distance
  G4double maxPathScale;	// Maximum ratio of trajectory vs. flight

  G4double flightAvg;		// Sum of flight distance for all steps
  G4double flightAvg2;		// Sum of squared flight distances (for RMS)
  G4double lastFlight;		// Last flight distance computed
  G4double lastRMS;		// Last flight distance RMS computed
  G4ThreeVector lastPos;	// Previous computed position

private:
  G4CMPTrackLimiter(const G4CMPTrackLimiter&);	// Copying is forbidden
  G4CMPTrackLimiter& operator=(const G4CMPTrackLimiter&);
};

#endif	/* G4CMPTrackLimiter_hh */
