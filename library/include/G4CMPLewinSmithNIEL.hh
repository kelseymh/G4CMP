/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLewinSmithNIEL.hh
/// \brief Non-ionizing energy loss calculation from Lewin & Smith 1996.
///
/// Computation of NIEL according to Lewin & Smith 1996, depending only
/// upon material properties.  Used by G4CMPEnergyPartition.
//
// $Id$
//
// 20190711  Michael Kelsey

#include "G4CMPVNIELPartition.hh"


class G4CMPLewinSmithNIEL : public G4CMPVNIELPartition {
public:
  G4LewinSmithNIEL() {;}
  virtual ~G4LewinSmithNIEL() {;}
  
  // return the fraction of the specified energy which will be deposited as NIEL
  // if an incoming particle with z1, a1 is stopped in the specified material
  // a1 is in atomic mass units, energy in native G4 energy units.
  //
  // NOTE: Projectile properties are ignored in Lewin & Smith calculation

  virtual G4double 
  PartitionNIEL(G4double energy, const G4Material *material, G4double Zin=0.,
		G4double Ain=0.) const;
};

#endif	/* G4CMPLewinSmithNIEL_hh */
