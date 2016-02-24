/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingDriftChamber.hh
/// \brief Definition of the ChannelingDriftChamber class
//
// $Id: a428d1182a10acb230bcab31560d634a3ea2f8f5 $
// --------------------------------------------------------------
//
#ifndef ChannelingDriftChamber_h
#define ChannelingDriftChamber_h 1

#include "G4VSensitiveDetector.hh"
#include "ChannelingDriftHit.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ChannelingDriftChamber : public G4VSensitiveDetector
{
  public:
      ChannelingDriftChamber(G4String name);
      virtual ~ChannelingDriftChamber();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      ChannelingDriftHitsCollection * fHitsCollection;
      G4int fHCID;
};

#endif

