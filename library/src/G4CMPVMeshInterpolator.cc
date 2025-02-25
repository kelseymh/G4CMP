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
//
// 20200914  Add function call to precompute potential gradients (field)
// 20240921  Add new Initialize() function to set tetra index cache
// 20250223  G4CMP-462: Avoid data race with worker thread Initialize()

#include "G4CMPVMeshInterpolator.hh"


// Replace values at mesh points without rebuilding tables

void G4CMPVMeshInterpolator::UseValues(const std::vector<G4double>& v) {
  if (!V.empty() && v.size() != V.size()) {
    G4cerr << "G4CMPVMeshInterpolator::UseValues ERROR Input vector v does"
	   << " not match existing mesh V." << G4endl;
    return;
  }

  V = v;
  FillGradients();	// Will call subclass implementation

#ifdef G4CMPTLI_DEBUG
  SavePoints(savePrefix+"_points.dat");
#endif
}


// Ensure that cached index is properly set

void G4CMPVMeshInterpolator::Initialize() {
  TetraIdx() = -1;
  if (TetraStart<0) TetraStart = FirstInteriorTetra();
  // FIXME: This may still cause a data race if a shared mesh instance
  //        is not instantiated by the master thread.
}
