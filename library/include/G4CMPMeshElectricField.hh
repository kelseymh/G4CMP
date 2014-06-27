// $Id$
//
// Wrapper class to process a numerically tabulated electric field mesh
// and use QHull to interpolate the potential and field at arbitrary
// points in the envelope of the mesh.  The input file format is fixed:
// each line consists of four floating-point values, x, y and z in meters,
// and voltage in volts.

#ifndef G4CMPMeshElectricField_h 
#define G4CMPMeshElectricField_h 1

#include "G4CMPTriLinearInterp.hh"
#include "G4ElectricField.hh"
#include <vector>

class G4CMPMeshElectricField : public G4ElectricField {
public:
  G4CMPMeshElectricField(const G4String& EpotFileName);
  //G4CMPMeshElectricField( G4double constEFieldVal );
  virtual ~G4CMPMeshElectricField() {;}

  virtual void GetFieldValue(const G4double Point[4], G4double *Efield) const;

  // Call through to interpolator (e.g., for use with FET code)
  virtual G4double GetPotential(const G4double Point[4]) const;

  // Copy constructor and assignment operator
  G4CMPMeshElectricField(const G4CMPMeshElectricField &p);
  G4CMPMeshElectricField& operator = (const G4CMPMeshElectricField &p);
  
private:
  G4CMPTriLinearInterp Interp;
  void BuildInterp(const G4String& EpotFileName);
};

namespace CDMS_Efield
{
  G4bool vector_comp(const std::vector<G4double>& p1,
		     const std::vector<G4double>& p2);
}

#endif	/* G4CMPMeshElectricField_h */
