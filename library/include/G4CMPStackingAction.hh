/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPStackingAction.hh
/// \brief Definition of the G4CMPStackingAction class
//
// $Id$
//
#ifndef G4CMPStackingAction_h
#define G4CMPStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4CMPProcessUtils.hh"

class G4Track;

class G4CMPStackingAction
  : public G4UserStackingAction, public G4CMPProcessUtils {
public:
  G4CMPStackingAction();
  virtual ~G4CMPStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);

protected:
  void SetPhononWaveVector(const G4Track* theTrack) const;
  void SetPhononVelocity(const G4Track* theTrack) const;

  void SetChargeCarrierValley(const G4Track* theTrack) const;
  void SetChargeCarrierMass(const G4Track* theTrack) const;
  void SetElectronEnergy(const G4Track* aTrack) const;
};

#endif /* G4CMPStackingAction_h */
