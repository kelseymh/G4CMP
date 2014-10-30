//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/include/G4CMPValleyTrackMap.hh
/// \brief Definition of the G4CMPValleyTrackMap base class
//
// $Id$
//
// 20131111  Move implementation of Clear() to .cc file
// 20141024  Support G4 before 10.0 by suppressing "G4ThreadLocal"

#ifndef G4CMPValleyTrackMap_h
#define G4CMPValleyTrackMap_h 1

#include "globals.hh"
#include <map>

#include "G4Version.hh"
#if G4VERSION_NUMBER < 1000
#define G4ThreadLocal
#endif

class G4Track;

class G4CMPValleyTrackMap {
public:
  static G4ThreadLocal G4CMPValleyTrackMap* theValleyTrackMap;

public:
  static G4CMPValleyTrackMap* GetValleyTrackMap();	// Synonyms for access
  static G4CMPValleyTrackMap* GetInstance() { return GetValleyTrackMap(); }

  // Update the wavevector for specified track, add track if not found
  void SetValley(const G4Track* track, G4int iv);
  void SetValley(const G4Track& track, G4int iv) { SetValley(&track, iv); }

  // Access current wavevector for specified track (NULL if doesn't exist)
  G4int GetValley(const G4Track* track) const;
  G4int GetValley(const G4Track& track) const { return GetValley(&track); }

  // Check if specified track is already loaded
  G4bool Find(const G4Track* track) const;
  G4bool Find(const G4Track& track) const { return Find(&track); }

  // Remove specified track from map (used by EndTracking)
  void RemoveTrack(const G4Track* track);

  void Clear();			// Remove all entries from map

private:
  typedef std::map<const G4Track*, G4int> TrkValleyMap;
  TrkValleyMap theMap;		// Associate tracks with vectors

private:
  G4CMPValleyTrackMap() { Clear(); }		// Ensure map is empty
  ~G4CMPValleyTrackMap() {;}
};
#endif	/* G4CMPValleyTrackMap_h */
