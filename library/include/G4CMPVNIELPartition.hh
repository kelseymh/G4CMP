/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPVNIELPartition.hh
/// \brief Definition of the G4CMPVNIELPartition base class
///
/// Abstract base class to define a "partition function" for non-ionizing
/// energy loss (NIEL) in material, which subclasses must implement.
//  Used by G4ScreenedNuclearRecoil and by G4CMPEnergyPartition.
//
// $Id$
//
// 20190711  Michael Kelsey

#ifndef G4CMPVNIELPartition_hh
#define G4CMPVNIELPartition_hh 1

#include "G4Types.hh"

class G4Material;


class G4CMPVNIELPartition {
public:
  G4VNIELPartition() {;}
  virtual ~G4VNIELPartition() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.
  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const = 0;
};

#endif	/* G4CMPVNIELPartition_hh */
