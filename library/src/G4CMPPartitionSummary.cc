/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPartitionSummary.hh
/// \brief Singleton container class to carry instances of G4CMPPartitionData
///   "hit" data, summarizing the output of G4CMPEnergyPartition.
///   
// $Id$
//
// 20200219  Michael Kelsey (TAMU) <kelsey@slac.stanford.edu>

#include "G4CMPPartitionSummary.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"


// Singleton construction

G4CMPPartitionSummary* G4CMPPartitionSummary::Instance() {
  static G4ThreadLocalSingleton<G4CMPPartitionSummary> theInstance;
  return theInstance.Instance();
}

// NOTE:  Data blocks use G4Allocator<>, which are deleted before TLSingletons

G4CMPPartitionSummary::~G4CMPPartitionSummary() {;}

// Delete previously collected data

void G4CMPPartitionSummary::clear() {
  for (auto& x: summary) { delete x; x=0; }
  summary.clear();
}

// Add new data record to collection, or restart if new event

void G4CMPPartitionSummary::insert(G4CMPPartitionData* data) {
  if (isNewEvent()) {			// Clean up previous event's data
    clear();
    currentEvent = getEventID();
  }

  if (data) summary.push_back(data);	// Only store valid data records
}

// Determine whether a new event has started

G4int G4CMPPartitionSummary::getEventID() {
  const G4Event* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  return (event ? event->GetEventID() : -1);
}

G4bool G4CMPPartitionSummary::isNewEvent() {
  return (getEventID() != currentEvent);
}
