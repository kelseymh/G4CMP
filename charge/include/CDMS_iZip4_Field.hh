#ifndef CDMS_iZip4_Field_h 
#define CDMS_iZip4_Field_h 1

#include "TriLinearInterp.hh"
#include "G4ElectricField.hh"

class CDMS_iZip4_Field : public G4ElectricField 
{
  public:
    CDMS_iZip4_Field( G4String EpotFileName ) : Interp(BuildInterp( EpotFileName )) {}
    //CDMS_iZip4_Field( G4double constEFieldVal );
    ~CDMS_iZip4_Field();

    void GetFieldValue(const G4double Point[4], G4double *Efield) const;
    void GetFieldValue(const G4double Point[3], G4double *Efield){
	    G4double NewPoint[4]={Point[0],Point[1],Point[2],0.0};
      this->GetFieldValue(NewPoint, Efield);
    };

    CDMS_iZip4_Field(const CDMS_iZip4_Field &p);
    CDMS_iZip4_Field& operator = (const CDMS_iZip4_Field &p);
    // Copy constructor and assignment operator
    
  private:
    TriLinearInterp Interp;
    TriLinearInterp BuildInterp( G4String EpotFileName );
};

namespace CDMS_Efield
{
  bool vector_comp( const vector<G4double>& p1, const vector<G4double>& p2 );
}


#endif
