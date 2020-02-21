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
