/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLindhardNIEL.hh
/// \brief Non-ionizing energy loss calculation from Lewin & Smith 1996.
///
/// Computation of NIEL according to Lindhard & Robinson, depending upon
/// both the material properties and the projectile mass and charge.
/// May be registered by client code into G4CMPEnergyPartition.
//
// $Id$
//
// 20190711  Michael Kelsey

#include "G4CMPVNIELPartition.hh"


class G4CMPLindhardNIEL : public G4CMPVNIELPartition {
public:
  G4LindhardNIEL() {;}
  virtual ~G4LindhardNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const;
};

#endif	/* G4CMPLindhardNIEL_hh */
