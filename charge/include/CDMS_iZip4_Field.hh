// $Id$

#ifndef CDMS_iZip4_Field_h 
#define CDMS_iZip4_Field_h 1

#include "G4CMPTriLinearInterp.hh"
#include "G4ElectricField.hh"
#include <vector>

class CDMS_iZip4_Field : public G4ElectricField {
public:
  CDMS_iZip4_Field(const G4String& EpotFileName);
  //CDMS_iZip4_Field( G4double constEFieldVal );
  virtual ~CDMS_iZip4_Field() {;}

  virtual void GetFieldValue(const G4double Point[4], G4double *Efield) const;
  virtual void GetFieldValue(const G4double Point[3], G4double *Efield);
  // NOTE:  This function is non-const ONLY because signatures are identical
  //        (the array dimension is not part of the function signature)

  // Copy constructor and assignment operator
  CDMS_iZip4_Field(const CDMS_iZip4_Field &p);
  CDMS_iZip4_Field& operator = (const CDMS_iZip4_Field &p);
  
private:
  G4CMPTriLinearInterp Interp;
  void BuildInterp(const G4String& EpotFileName);
};

namespace CDMS_Efield
{
  G4bool vector_comp(const std::vector<G4double>& p1,
		     const std::vector<G4double>& p2);
}

#endif	/* CDMS_iZip4_Field_h */
