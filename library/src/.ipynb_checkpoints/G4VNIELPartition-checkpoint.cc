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
// $Id: 875173f1491afd2935ee7cdeb290662c0ea5e837 $
//
// 20191211  Michael Kelsey

#include "G4VNIELPartition.hh"
#include "G4Element.hh"
#include "G4Material.hh"


// Computed weighted average of elemental constituents of material

G4double G4VNIELPartition::GetEffectiveZ(const G4Material *material) const {
  if (!material) {
    G4Exception ("G4VNIELPartition::GetEffectiveZ()", "G4CMP1001",
		 FatalException, "Called with null material pointer");
    return 0.;
  }

  // Accumulate weighted Z for all elements in material
  const G4Element* element = 0;
  G4double zEff = 0.;
  for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
    element = material->GetElement(i);
    zEff += element->GetZ() * material->GetVecNbOfAtomsPerVolume()[i];
  }

  // Normalize to total material
  zEff /= material->GetTotNbOfAtomsPerVolume();

  return zEff;
}

G4double G4VNIELPartition::GetEffectiveA(const G4Material *material) const {
  if (!material) {
    G4Exception ("G4VNIELPartition::GetEffectiveA()", "G4CMP1002",
		 FatalException, "Called with null material pointer");
    return 0.;
  }

  // Accumulate weighted A for all elements in material
  const G4Element* element = 0;
  G4double aEff = 0.;
  for (unsigned int i = 0; i < material->GetNumberOfElements(); ++i) {
    element = material->GetElement(i);
    aEff += element->GetA() * material->GetVecNbOfAtomsPerVolume()[i];
  }

  // Normalize to total material
  aEff /= material->GetTotNbOfAtomsPerVolume();

  return aEff;
}
