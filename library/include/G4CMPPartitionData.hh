/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPartitionData.hh
/// \brief Definition of the G4CMPPartitionData "hit" (G4VHit) container
///   Simple container class used by G4CMPEnergyPartition to record
///   some of the internal quantities used to generate an event.  Stored
///   into the G4Event via G4CMPPartitionCollection (G4VHitsCollection) 
///   for access by user applications.
///   
// $Id$
//
// 20200217  Michael Kelsey (TAMU) <kelsey@slac.stanford.edu>

#ifndef G4CMPPartitionData_hh
#define G4CMPPartitionData_hh 1

#include "G4Types.hh"
#include "G4Allocator.hh"
#include "G4VHit.hh"


class G4CMPPartitionData : public G4VHit {
public:
  G4CMPPartitionData() : G4VHit() {;}
  virtual ~G4CMPPartitionData() {;}

  // Memory allocation using G4Allocator<>; implemented below
  inline void* operator new(size_t);
  inline void  operator delete(void*);

public:		// Simple container, provide direct access to information
  G4double trueEnergy;		// Input total energy deposited (dE/dx)
  G4double trueNIEL;		// Input non-ionizing energy deposited
  G4double lindhardYield;	// Computed value, or 1 - NIEL/Energy
  G4double phononEnergy;	// Energy assigned for primary phonons
  G4double phononGenerated;	// Weighed sum of generated primary phonons
  G4double chargeEnergy;	// Energy assigned to charge carriers
  G4double chargeGenerated;	// Fano-fluctuated energy in charge carriers
  G4int numberOfPairs;		// Number of e/h pairs created (downsampled)

public:
  G4CMPPartitionData(const G4CMPPartitionData&) = default;
  G4CMPPartitionData(G4CMPPartitionData&& = default;
  G4CMPPartitionData& operator=(const G4CMPPartitionData&) = default;
  G4CMPPartitionData& operator=(G4CMPPartitionData&& = default;
};


// Data and memory management

#if defined G4DIGI_ALLOC_EXPORT
  extern G4DLLEXPORT 
  G4ThreadLocal G4Allocator<G4CMPPartitionData>* G4CMPPartitionData_Allocator;
#else
  extern G4DLLIMPORT 
  G4ThreadLocal G4Allocator<G4CMPPartitionData>* G4CMPPartitionData_Allocator;
#endif

// Allocation operators should be inlined for efficiency

inline void* G4CMPPartitionData::operator new(size_t) {
  if (!G4CMPPartitionData_Allocator)
    G4CMPPartitionData_Allocator = new G4Allocator<G4CMPPartitionData>;
  return (void*) G4CMPPartitionData_Allocator->MallocSingle();
}

inline void G4CMPPartitionData::operator delete(void* aHit) {
  G4CMPPartitionData_Allocator->FreeSingle((G4CMPPartitionData*) aHit);
}

#endif	/* G4CMPPartitionData_hh */
