/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChannelingTrackingAction_h
#define ChannelingTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class ChannelingTrackingAction : public G4UserTrackingAction {
public:
  void PreUserChannelingTrackingAction(const G4Track*);
  void PostUserChannelingTrackingAction(const G4Track*);
};

#endif










