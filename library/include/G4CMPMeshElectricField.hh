/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// Wrapper class to process a numerically tabulated electric field mesh
// and use Qhull to interpolate the potential and field at arbitrary
// points in the envelope of the mesh.  The input file format is fixed:
// each line consists of four floating-point values, x, y and z in meters,
// and voltage in volts.
//
// 20150122  Move vector_comp into class definition as static function.
// 20170823  Add scaling factor as optional constructor argument.
// 20180904  Add constructor to take precreated mesh and tetrahedra.
// 20180924  TriLinearInterp should be a pointer, to break dependency.
// 20190226  Provide access to TriLinearInterp object, and ctor assignment

#ifndef G4CMPMeshElectricField_h 
#define G4CMPMeshElectricField_h 1

#include "G4ElectricField.hh"
#include "G4ThreeVector.hh"
#include <array>
#include <vector>

class G4CMPTriLinearInterp;


class G4CMPMeshElectricField : public G4ElectricField {
public:
  G4CMPMeshElectricField(const G4String& EPotFileName, G4double Vscale=1.);

  G4CMPMeshElectricField(const std::vector<std::array<G4double,3> >& xyz,
			 const std::vector<G4double>& v,
			 const std::vector<std::array<G4int,4> >& tetra);

  G4CMPMeshElectricField(const G4CMPTriLinearInterp& tli);

  // Copy constructor and assignment operator
  G4CMPMeshElectricField(const G4CMPMeshElectricField &p);
  G4CMPMeshElectricField& operator=(const G4CMPMeshElectricField &p);

  virtual ~G4CMPMeshElectricField();

  // Returns field vector(s) at location, required by base class
  virtual void GetFieldValue(const G4double Point[3], G4double *Efield) const;

  // Call through to interpolator (e.g., for use with FET code)
  virtual G4double GetPotential(const G4double Point[3]) const;

  // Get access to mesh interpolator for client access or copying
  const G4CMPTriLinearInterp* GetInterpolator() const { return Interp; }

  // Sorting operator (compares x, y, z in sequence)
  static G4bool vector_comp(const std::array<G4double, 4>& p1,
			    const std::array<G4double, 4>& p2);

private:
  G4CMPTriLinearInterp* Interp;

  void BuildInterp(const G4String& EPotFileName, G4double Vscale=1.);

  void BuildInterp(const std::vector<std::array<G4double,3> >& xyz,
		   const std::vector<G4double>& v,
		   const std::vector<std::array<G4int,4> >& tetra);
};

#endif	/* G4CMPMeshElectricField_h */
