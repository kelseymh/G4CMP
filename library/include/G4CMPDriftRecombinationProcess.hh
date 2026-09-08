/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPDriftRecombinationProcess.hh
/// \brief Class definition for charge (electron-hole) recombination
//
// 20170802  Add EnergyPartition to handle phonon production
// 20260514  Add line breaks at 80 columns, for readability
// 20260514  G4CMP-517 -- Support recombination in flight with energy gain
//	       calculations.

#ifndef G4CMPDriftRecombinationProcess_h
#define G4CMPDriftRecombinationProcess_h 1

#include "G4CMPVDriftProcess.hh"
#include "globals.hh"

class G4CMPEnergyPartition;


class G4CMPDriftRecombinationProcess : public G4CMPVDriftProcess {
public:
  G4CMPDriftRecombinationProcess(const G4String& name = "G4CMPChargeRecombine",
                                 G4CMPProcessSubType type = fChargeRecombine);
  virtual ~G4CMPDriftRecombinationProcess();

  // No copying/moving
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&) = delete;
  G4CMPDriftRecombinationProcess(G4CMPDriftRecombinationProcess&&) = delete;

  G4CMPDriftRecombinationProcess&
  operator=(const G4CMPDriftRecombinationProcess&) = delete;

  G4CMPDriftRecombinationProcess&
  operator=(const G4CMPDriftRecombinationProcess&&) = delete;

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&)
    override;

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
    override;

  // Flag if track is eligible for recombination: stopped or below threshold
  G4bool ReadyToRecombine(const G4Track& aTrack) const;

  // Decide if track gains enough energy for NTL emission before surface
  // NOTE: Will include "turn around" if electric field is in use
  G4bool EnergyGainToSurface(const G4Track& aTrack) const;
  G4bool LukeBeforeSurface(const G4Track& aTrack) const;

  // Compute distance and direction to surface, using field acceleration
  G4ThreeVector VectorToSurface(const G4Track& aTrack) const;

  // Compute acceleration direction from E-field, accounting for valleys
  G4ThreeVector GetAcceleration(const G4ThreeVector& Efield) const;

private:
  G4CMPEnergyPartition* partitioner;
};

#endif	/* G4CMPDriftRecombinationProcess_h */
