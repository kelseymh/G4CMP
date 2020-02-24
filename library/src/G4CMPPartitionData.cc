/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPPartitionData.cc
/// \brief Definition of the G4CMPPartitionData "hit" (G4VHit) container
///   Simple container class used by G4CMPEnergyPartition to record
///   some of the internal quantities used to generate an event.  Stored
///   into the G4Event via G4CMPPartitionCollection (G4VHitsCollection) 
///   for access by user applications.
///   
// $Id$
//
// 20200218  Michael Kelsey (TAMU) <kelsey@slac.stanford.edu>

#include "globals.hh"
#include "G4CMPPartitionData.hh"
#include "G4SystemOfUnits.hh"


// Initialize allocator

G4ThreadLocal G4Allocator<G4CMPPartitionData>* G4CMPPartitionData_Allocator=0;

// Constructor must initialize everything to zero

G4CMPPartitionData::G4CMPPartitionData()
  : G4VHit(), totalEnergy(0.), truedEdx(0.), trueNIEL(0.), lindhardYield(0.),
    FanoFactor(0.), chargeEnergy(0.), chargeGenerated(0.), numberOfPairs(0),
    phononEnergy(0.), phononGenerated(0.), numberOfPhonons(0) {;}

// Report contents for diagnostic purposes

void G4CMPPartitionData::Print() {	// FIXME: Base class is non-const
  G4cout << "G4CMPPartitionData:"
	 << "\n Total energy deposit " << totalEnergy/eV << " eV"
	 << "\n Ionization energy (dE/dx) " << truedEdx/eV << " eV"
	 << "\n Non-ionizing energy " << trueNIEL/eV << " eV"
	 << "\n Lindhard yield " << lindhardYield
	 << "\n True charge energy " << chargeEnergy/eV << " eV"
	 << "\n Fano factor " << FanoFactor
	 << "\n Generated charge energy " << chargeGenerated/eV << " eV"
	 << "\n Number of charge pairs " << numberOfPairs
	 << "\n True phonon energy " << phononEnergy/eV << " eV"
	 << "\n Generated phonon energy " << phononGenerated/eV << " eV"
	 << "\n Number of phonons " << numberOfPhonons
	 << G4endl;
}
