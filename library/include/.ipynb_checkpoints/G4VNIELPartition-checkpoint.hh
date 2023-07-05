/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4VNIELPartition.hh
/// \brief Definition of the G4VNIELPartition base class
///
/// Abstract base class to define a "partition function" for non-ionizing
/// energy loss (NIEL) in material, which subclasses must implement.
//  Used by G4ScreenedNuclearRecoil and by G4CMPEnergyPartition.
//
// $Id: 14f4e2c177535b52297b78008bc7ab103d6b1885 $
//
// 20190711  Michael Kelsey
// 20191211  Add functions to compute effective Z and A of composite material

#ifndef G4VNIELPartition_hh
#define G4VNIELPartition_hh 1

#include "G4Types.hh"

class G4Material;


class G4VNIELPartition {
public:
  G4VNIELPartition() {;}
  virtual ~G4VNIELPartition() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.
  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const = 0;

protected:
  G4double GetEffectiveZ(const G4Material *material) const;
  G4double GetEffectiveA(const G4Material *material) const;
};

#endif	/* G4CMPVNIELPartition_hh */
