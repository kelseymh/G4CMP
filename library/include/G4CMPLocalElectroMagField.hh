// $Id$
//
// Wrapper class for G4ElectroMagneticField objects which handles transforming
// between global and local coordinates for input and output queries.

#ifndef G4CMPLocalElectroMagField_hh
#define G4CMPLocalElectroMagField_hh 1

#include "G4ElectroMagneticField.hh"
#include "G4AffineTransform.hh"


class G4CMPLocalElectroMagField : public G4ElectroMagneticField {
public:
  G4CMPLocalElectroMagField(const G4ElectroMagneticField* theField)
    : G4ElectroMagneticField(), localField(theField) {;}

  virtual ~G4CMPLocalElectroMagField() {;}

  G4CMPLocalElectroMagField(const G4CMPLocalElectroMagField& rhs)
    : G4ElectroMagneticField(rhs), localField(rhs.localField) {;}

  G4CMPLocalElectroMagField& operator=(const G4CMPLocalElectroMagField& rhs) {
    G4ElectroMagneticField::operator=(rhs);
    localField = rhs.localField;
    return *this;
  }

  // Specify local-to-global transformation before field call
  // Typically, this will be called from FieldManager::ConfigureForTrack()
  virtual void SetTransforms(const G4AffineTransform& lToG);

  // This function transforms Point[0..2] to local coordinates on input,
  // and Field[0..2], Field[3..5] from local to global coordinate on output
  virtual void GetFieldValue(const G4double Point[4], G4double *BEfield) const;

protected:
  void GetLocalPoint(const G4double Point[4]) const;
  void CopyLocalToGlobalVector(G4int index, G4double* gbl) const;

private:
  const G4ElectroMagneticField* localField;	// Local-coord implementation
  G4AffineTransform fLocalToGlobal;
  G4AffineTransform fGlobalToLocal;

  mutable G4ThreeVector vec;		// Internal buffers to reduce memory
  mutable G4double localP[4];
  mutable G4double localF[6];
};

#endif /* G4CMPLocalElectroMagField_hh */
