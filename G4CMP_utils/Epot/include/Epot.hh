#ifndef DetectorField_h
#define DetectorField_h
#include <vector>
#include "globals.hh"
#include "G4ElectricField.hh"
#include "Interpolation3D.hh"

using std::vector;

class DetectorField : public G4ElectricField 
{
    public:
        DetectorField(G4String filename);
    	~DetectorField();
        void GetFieldValue(const G4double Point[4], G4double *Efield) const;
        void GetFieldValue(const G4double Point[3], G4double *Efield){
		    G4double NewPoint[4]={Point[0],Point[1],Point[2],0.0};
		    this->GetFieldValue(NewPoint, Efield);
	    };
        inline void GetX(vector<vector<G4double> >& mat){mat=xyzMat;}
        inline void GetV(vector<G4double>& volt){volt=V;}

    private:
        vector<vector<G4double> > xyzMat;
	    vector<G4double> V;
	    G4double xmin,xmax,ymin,ymax,zmin,zmax;
        Interpolation3D Interpolation;
};
#endif
