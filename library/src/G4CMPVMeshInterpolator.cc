/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// G4CMPVMeshInterpolator:  Base interface class for BiLinear and TriLinear
// interpolators, used by G4CMPMeshElectricField.  This base class allows
// use of a single pointer with a type selected at runtime for the field
// model.  The actual interpolation functionality is implemented entirely
// in the concrete subclasses.

#include "G4CMPVMeshInterpolator.hh"


// Replace values at mesh points without rebuilding tables

void G4CMPVMeshInterpolator::UseValues(const std::vector<G4double>& v) {
  if (!V.empty() && v.size() != V.size()) {
    G4cerr << "G4CMPVMeshInterpolator::UseValues ERROR Input vector v does"
	   << " not match existing mesh V." << G4endl;
    return;
  }

  staleCache = true;
  V = v;

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
#endif
}

