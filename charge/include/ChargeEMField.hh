#ifndef ChargeEMField_h
#define ChargeEMField_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class ChargeFieldMessenger;

class ChargeEMField
{
public:
  ChargeEMField(G4LogicalVolume* logVol);
  ~ChargeEMField();
  void Build();

  void SetEPotFileName(G4String filename) {EPotFile = filename;}
  void UseConstantField(G4bool constField) {UseConstEField = constField;}
  void SetFieldMagnitude(G4double mag) {ConstEFieldMag = mag;}
  void SetFieldDirection(G4ThreeVector direct) {ConstEFieldDir = direct.unit();}

  private:
    ChargeFieldMessenger* Messenger;
    G4String EPotFile;
    G4bool UseConstEField;
    G4double ConstEFieldMag;
    G4ThreeVector ConstEFieldDir;
    G4LogicalVolume* DetectorVolume;
};

#endif
