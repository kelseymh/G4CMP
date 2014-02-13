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
/// \file library/src/G4CMPValleyTrackMap.hh
/// \brief Implementation of the G4CMPValleyTrackMap base class
//
// $Id$

#include "G4CMPValleyTrackMap.hh"
#include "G4Track.hh"
#include <map>


// Singleton instance, must be created at first request

G4ThreadLocal G4CMPValleyTrackMap* G4CMPValleyTrackMap::theValleyTrackMap = 0;


// Pseudo-constructor creates singleton instance

G4CMPValleyTrackMap* G4CMPValleyTrackMap::GetValleyTrackMap() {
  if (!theValleyTrackMap) theValleyTrackMap = new G4CMPValleyTrackMap;
  return theValleyTrackMap;
}

void G4CMPValleyTrackMap::Clear() {
  theMap.clear();			// Remove all entries from map
}



// Check if specified track is already loaded

G4bool G4CMPValleyTrackMap::Find(const G4Track* track) const {
  return (!track || theMap.find(track) != theMap.end());
}


// Remove specified track from map (used by EndTracking)

void G4CMPValleyTrackMap::RemoveTrack(const G4Track* track) {
  TrkValleyMap::iterator entry = theMap.find(track);
  if (entry != theMap.end()) theMap.erase(entry);
}


// Update the wavevector for specified track, add track if non-existent

void G4CMPValleyTrackMap::SetValley(const G4Track* track, G4int iv) {
  if (track) theMap[track] = iv;
}


// Access current wavevector for specified track (NULL if doesn't exist)

G4int G4CMPValleyTrackMap::GetValley(const G4Track* track) const {
  TrkValleyMap::const_iterator entry = theMap.find(track);
  return (entry != theMap.end()) ? entry->second : -1;
}
  

