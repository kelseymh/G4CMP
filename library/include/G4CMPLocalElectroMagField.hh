/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// Wrapper class for G4ElectroMagneticField objects which handles transforming
// between global and local coordinates for input and output queries.
//
// 20180525  Provide accessor to underlying (local coordinate) field.
// 20180711  Store local geometry associated with field, and accessor to
//		optionally interpolate potential
// 20210902  Add verbosity flag set in constructing client code

#ifndef G4CMPLocalElectroMagField_hh
#define G4CMPLocalElectroMagField_hh 1

#include "G4ElectroMagneticField.hh"
#include "G4AffineTransform.hh"

class G4VSolid;


class G4CMPLocalElectroMagField : public G4ElectroMagneticField {
public:
  G4CMPLocalElectroMagField(const G4ElectroMagneticField* theField)
    : G4ElectroMagneticField(), localField(theField), theSolid(0),
      verboseLevel(0) {;}

  virtual ~G4CMPLocalElectroMagField() {;}

  G4CMPLocalElectroMagField(const G4CMPLocalElectroMagField& rhs)
    : G4ElectroMagneticField(rhs), localField(rhs.localField),
      theSolid(rhs.theSolid), fLocalToGlobal(rhs.fLocalToGlobal),
      fGlobalToLocal(rhs.fGlobalToLocal), verboseLevel(rhs.verboseLevel) {;}

  G4CMPLocalElectroMagField& operator=(const G4CMPLocalElectroMagField& rhs) {
    G4ElectroMagneticField::operator=(rhs);
    localField = rhs.localField;
    theSolid = rhs.theSolid;
    fLocalToGlobal = rhs.fLocalToGlobal;
    fGlobalToLocal = rhs.fGlobalToLocal;
    verboseLevel = rhs.verboseLevel;

    return *this;
  }

  // Turn on diagnostic output
  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
  G4int GetverboseLevel() const { return verboseLevel; }

  // Specify volume and local-to-global transformation before field call
  // Typically, this will be called from FieldManager::ConfigureForTrack()
  virtual void SetTransforms(const G4AffineTransform& lToG);
  virtual void SetGeometry(const G4VSolid* shape) { theSolid = shape; }

  // Call-through to user-defined field, since properties aren't changed
  virtual G4bool DoesFieldChangeEnergy() const {
    return localField->DoesFieldChangeEnergy();
  }

  // This function transforms Point[0..2] to local coordinates on input,
  // and Field[0..2], Field[3..5] from local to global coordinate on output
  virtual void GetFieldValue(const G4double Point[4], G4double *BEfield) const;

  // This function takes a GLOBAL position and returns the potential at
  // that point, using the stored volume to integrate up to the surfaces.
  virtual G4double GetPotential(const G4double Point[4]) const;

  // Provide client code with access to original local-coordinate field
  // This is mainly useful for field subclasses with extended interfaces
  const G4ElectroMagneticField* GetLocalField() const { return localField; }

protected:
  void GetLocalPoint(const G4double Point[4]) const;
  void CopyLocalToGlobalVector(G4int index, G4double* gbl) const;

private:
  const G4ElectroMagneticField* localField;	// Local-coord implementation
  const G4VSolid* theSolid;		// Shape to which field is attached
  G4AffineTransform fLocalToGlobal;
  G4AffineTransform fGlobalToLocal;

  G4int verboseLevel;			// For diagnostic output

  mutable G4ThreeVector vec;		// Internal buffers to reduce memory
  mutable G4double localP[4];
  mutable G4double localF[6];
};

#endif /* G4CMPLocalElectroMagField_hh */
