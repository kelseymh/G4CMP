#ifndef ChargeEMField_h
#define ChargeEMField_h 1

class G4LogicalVolume;

class ChargeEMField
{
public:
  ChargeEMField(G4LogicalVolume* logVol);
  ~ChargeEMField();
};

#endif
