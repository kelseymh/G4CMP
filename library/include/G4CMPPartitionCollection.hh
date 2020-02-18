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

#ifndef G4CMPPartitionCollection_hh
#define G4CMPPartitionCollection_hh 1

#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4CMPPartitionData.hh"


// Provide convenient name for hits container within collection
typedef std::vector<G4CMPPartitionData*> G4CMPPartitionVector;

class G4CMPPartitionCollection : public G4THitsCollection<G4CMPPartitionData> {
public:		// Provide references to hardcoded names for client code
  static const char* sdName;	// "G4CMP"
  static const char* colName;	// "G4CMPPartition"
  static const char* keyName;	// "G4CMP/G4CMPPartition"

public:
  G4CMPPartitionCollection()
    : G4THitsCollection<G4CMPPartitionData>(sdName,colName) {;}
  virtual ~G4CMPPartitionCollection() {;}

  // Memory allocation using G4Allocator<>; implemented below
  inline void* operator new(size_t);
  inline void  operator delete(void*);
};


// Data and memory management

#if defined G4DIGI_ALLOC_EXPORT
  extern G4DLLEXPORT G4ThreadLocal G4Allocator<G4CMPPartitionCollection>* G4CMPPartColl_Allocator;
#else
  extern G4DLLIMPORT G4ThreadLocal G4Allocator<G4CMPPartitionCollection>* G4CMPPartColl_Allocator;
#endif

// Allocation operators should be inlined for efficiency

inline void* G4CMPPartitionCollection::operator new(size_t) {
  if (!G4CMPPartColl_Allocator)
    G4CMPPartColl_Allocator = new G4Allocator<G4CMPPartitionCollection>;
  return (void*) G4CMPPartColl_Allocator->MallocSingle();
}

inline void G4CMPPartitionCollection::operator delete(void* aHit) {
  G4CMPPartColl_Allocator->FreeSingle((G4CMPPartitionCollection*) aHit);
}

#endif	/* G4CMPPartitionCollection_hh */
