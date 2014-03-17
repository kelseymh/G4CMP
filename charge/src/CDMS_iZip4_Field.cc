#include "CDMS_iZip4_Field.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <fstream>

using std::vector;


CDMS_iZip4_Field::CDMS_iZip4_Field(const G4String& EpotFileName)
  : G4ElectricField() {
  BuildInterp(EpotFileName);
}

void CDMS_iZip4_Field::BuildInterp(const G4String& EpotFileName) {
  G4cout << "CDMS_iZip4_Field::Constructor: Creating Electric Field" << G4endl;
  vector<vector<G4double> > tempX;

  vector<G4double> temp = vector<G4double>(4, 0);
  G4double x,y,z,v;

  std::ifstream epotFile(EpotFileName, std::ios::in);
  while (!epotFile.eof())
  {
    epotFile >> x >> y >> z >> v;
    temp[0] = x*m;
    temp[1] = y*m;
    temp[2] = z*m - 1.27*cm;	// *** FIXME:  THIS IS THICKNESS OF CDMS iZIP!
    temp[3] = v*volt;
    tempX.push_back(temp);
  }
  epotFile.close();

  std::sort(tempX.begin(),tempX.end(), CDMS_Efield::vector_comp);

  vector<vector<G4double> > X(tempX.size(),vector<G4double>(3,0));
  vector<G4double> V(tempX.size(),0);

  for (size_t ii = 0; ii < tempX.size(); ++ii)
  {
    X[ii][0] = tempX[ii][0];
    X[ii][1] = tempX[ii][1];
    X[ii][2] = tempX[ii][2];
    V[ii] = tempX[ii][3];
  }

  Interp.UseMesh(X, V);
}


G4bool CDMS_Efield::vector_comp(const vector<G4double>& p1,
				const vector<G4double>& p2) {
  if (p1[0] < p2[0])
    return true;
  else if (p2[0] < p1[0])
    return false;
  else if (p1[1] < p2[1])
    return true;
  else if (p2[1] < p1[1])
    return false;
  else if (p1[2] < p2[2])
    return true;
  else if (p2[2] < p1[2])
    return false;
  else
    return false;
}

void CDMS_iZip4_Field::GetFieldValue(const G4double Point[3],
				     G4double *Efield) {
  G4double tempPoint[4] = {Point[0], Point[1], Point[2], 0.0};
  this->GetFieldValue(tempPoint, Efield);
};

void CDMS_iZip4_Field::GetFieldValue(const G4double Point[4],
				     G4double *Efield) const {
  Interp.GetField(Point,Efield);
}

CDMS_iZip4_Field::CDMS_iZip4_Field(const CDMS_iZip4_Field &p)
  : G4ElectricField(p), Interp(p.Interp) {;}

CDMS_iZip4_Field& CDMS_iZip4_Field::operator=(const CDMS_iZip4_Field &p) {
  if (this != &p) {				// Only copy if not self
    G4ElectricField::operator=(p);		// Call through to base
    Interp = p.Interp;
  }

  return *this;
}
