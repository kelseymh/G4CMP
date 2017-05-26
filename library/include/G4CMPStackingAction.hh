/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPStackingAction.hh
/// \brief Definition of the G4CMPStackingAction class
//
// $Id$
//
// 20170525  M. Kelsey -- Add default "rule of five" copy/move operators

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
  void SetPhononVelocity(const G4Track* theTrack) const;

  void SetChargeCarrierMass(const G4Track* theTrack) const;
  void SetElectronEnergy(const G4Track* aTrack) const;

public:
  G4CMPStackingAction(const G4CMPStackingAction&) = default;
  G4CMPStackingAction(G4CMPStackingAction&&) = default;
  G4CMPStackingAction& operator=(const G4CMPStackingAction&) = default;
  G4CMPStackingAction& operator=(G4CMPStackingAction&&) = default;

};

#endif /* G4CMPStackingAction_h */
