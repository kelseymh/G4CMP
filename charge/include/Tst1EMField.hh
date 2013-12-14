#ifndef Tst1EMField_h
#define Tst1EMField_h 1

#include "EqEMFieldXtal.hh"
#include "G4UniformElectricField.hh"

class G4LogicalVolume;

class Tst1EMField
{
public:
  Tst1EMField(G4LogicalVolume* logVol);
  ~Tst1EMField();
};


#endif
