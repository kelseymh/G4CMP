/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPartitionCollection.hh
/// \brief Definition of the G4VHitsCollection subclass to carry instances
///   of G4CMPPartitionData "hit" data.
///   
// $Id$
//
// 20200217  Michael Kelsey (TAMU) <kelsey@slac.stanford.edu>

#include "G4CMPPartitionCollection.hh"


// Initialize hardcoded names for collection

const char* G4CMPPartitionCollection::sdName  = "G4CMP";
const char* G4CMPPartitionCollection::colName = "G4CMPPartition";
const char* G4CMPPartitionCollection::keyName = "G4CMP/G4CMPPartition";

// Initialize allocator

G4ThreadLocal G4Allocator<G4CMPPartitionCollection>* G4CMPPartColl_Allocator=0;
