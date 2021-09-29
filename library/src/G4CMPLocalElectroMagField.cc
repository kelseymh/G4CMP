/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// Wrapper class for G4ElectroMagneticField objects which handles transforming
// between global and local coordinates for input and output queries.
//
// 20180711  Provide interpolator to return potential at point in volume,
//	       assuming "mid-plane" is at ground.

#include "G4CMPLocalElectroMagField.hh"
#include "G4CMPMeshElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4VSolid.hh"
#include <algorithm>

// Specify local-to-global transformation before field call
// Typically, this will be called from FieldManager::ConfigureForTrack()

void G4CMPLocalElectroMagField::SetTransforms(const G4AffineTransform& lToG) {
  if (verboseLevel>2) {
    G4cout << "LocalEMField::SetTransforms ltoG "
	   << " trans " << lToG.NetTranslation() << " rot "
	   << lToG.NetRotation().delta()/CLHEP::deg << " deg about "
	   << lToG.NetRotation().axis() << G4endl;
  }

  fGlobalToLocal = fLocalToGlobal = lToG;
  fGlobalToLocal.Invert();

  if (verboseLevel>2) {
    G4ThreeVector glob111 = G4ThreeVector(1,1,1);
    G4ThreeVector loc111 = fGlobalToLocal.TransformPoint(glob111);
    G4cout << " gToL moves " << glob111 << " -> " << loc111 << G4endl
	   << " lToG moves " << loc111*10. << " -> "
	   << fLocalToGlobal.TransformPoint(loc111*10.) << G4endl;
  }
}


// This function transforms Point[0..2] to local coordinates on input,
// and Field[3..5] from local to global coordinate on output
void G4CMPLocalElectroMagField::GetFieldValue(const G4double Point[4],
					      G4double *BEfield) const {
  GetLocalPoint(Point);

  // Get local field vector(s) using local position
  std::fill(localF, localF+6, 0.);
  localField->GetFieldValue(localP, localF);

  // Transform local magnetic and electric fields to global
  CopyLocalToGlobalVector(0, BEfield);
  CopyLocalToGlobalVector(3, BEfield);
}


// Copy input point to vector and rotate to local coordinates

void G4CMPLocalElectroMagField::GetLocalPoint(const G4double Point[4]) const {
  vec.set(Point[0], Point[1], Point[2]);
  fGlobalToLocal.ApplyPointTransform(vec);

  if (verboseLevel>2) {
    G4cout << "LocalEMField::GetLocalPoint (" << Point[0]
	   << "," << Point[1] << "," << Point[2] << ") -> " << vec << G4endl;
  }

  // Fill local point buffer and initialize field
  localP[0] = vec.x();
  localP[1] = vec.y();
  localP[2] = vec.z();
  localP[3] = Point[3];		// Time coordinate
}

// Transform local vector (array) to global coordinate system

void G4CMPLocalElectroMagField::
CopyLocalToGlobalVector(G4int index, G4double* gbl) const {
  vec.set(localF[index+0], localF[index+1], localF[index+2]);
  fLocalToGlobal.ApplyAxisTransform(vec);

  if (verboseLevel>2) {
    G4cout << "LocalEMField::CopyLocalToGlobalVector " << index
	   << " (" << localF[index+0] << "," << localF[index+1] << ","
	   << localF[index+2] << ") -> " << vec << G4endl;
  }

  gbl[index+0] = vec.x();
  gbl[index+1] = vec.y();
  gbl[index+2] = vec.z();
}


// This function takes a GLOBAL position and returns the potential at
// that point, using the stored volume to integrate up to the surfaces.

G4double G4CMPLocalElectroMagField::
GetPotential(const G4double Point[4]) const {
  GetLocalPoint(Point);

  // Mesh field can return potential directly; no work needed
  const G4CMPMeshElectricField* mfield =
    dynamic_cast<const G4CMPMeshElectricField*>(localField);
  if (mfield) {
    if (verboseLevel>2) {
      G4cout << "LocalEMField::GetPotential mesh returns " 
	     << mfield->GetPotential(localP) << G4endl;
    }

    return mfield->GetPotential(localP);
  }

  if (!theSolid) return 0.;		// Can't integrate without boundaries

  // Uniform field doesn't require stepwise integration
  const G4UniformElectricField* ufield =
    dynamic_cast<const G4UniformElectricField*>(localField);
  if (ufield) {
    ufield->GetFieldValue(localP, localF);
    G4ThreeVector evec(localF[3],localF[4],localF[5]);
    G4ThreeVector pos(localP[0],localP[1],localP[2]);

    G4ThreeVector e0 = evec.unit();
    G4double toVpos = theSolid->DistanceToOut(pos, -e0);
    G4double toVneg = theSolid->DistanceToOut(pos, e0);

    if (verboseLevel>2) {
      G4cout << "LocalEMField::GetPotential pos " << pos << " e0 " << e0
	     << " toVpos " << toVpos << " toVneg " << toVneg
	     << " emag " << evec.mag()
	     << " : V = " << 0.5*(toVpos-toVneg)*evec.mag()
	     << G4endl;
    }

    return 0.5*(toVpos-toVneg)*evec.mag();	// [-V/2,V/2] interpolation
  }

  // Arbitrary field configurations must be integrated
  return 0.;
}
