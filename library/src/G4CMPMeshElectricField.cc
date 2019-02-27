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
// 20150122  Use verboseLevel instead of compiler flag for debugging; move
//	     vector_comp into class (static).
// 20170823  Add scaling factor as optional constructor argument.
// 20180525  Use new "quiet" argument to suppress "outside of hull" warnings
// 20180904  Add constructor to take precreated mesh and tetrahedra.
// 20180924  TriLinearInterp should be a pointer, to break dependency.
// 20190226  Provide access to TriLinearInterp object, and ctor assignment

#include "G4CMPMeshElectricField.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPTriLinearInterp.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

using std::array;
using std::vector;


// Constructors

G4CMPMeshElectricField::
G4CMPMeshElectricField(const G4String& EPotFileName, G4double Vscale)
  : G4ElectricField(), Interp(new G4CMPTriLinearInterp) {
  BuildInterp(EPotFileName, Vscale);
}

G4CMPMeshElectricField::
G4CMPMeshElectricField(const vector<array<G4double,3> >& xyz,
		       const vector<G4double>& v,
		       const vector<array<G4int,4> >& tetra)
  : G4ElectricField(), Interp(new G4CMPTriLinearInterp) {
  BuildInterp(xyz, v, tetra);
}

G4CMPMeshElectricField::G4CMPMeshElectricField(const G4CMPMeshElectricField &p)
  : G4ElectricField(p), Interp(new G4CMPTriLinearInterp(*p.Interp)) {;}

G4CMPMeshElectricField::G4CMPMeshElectricField(const G4CMPTriLinearInterp& tli)
  : G4ElectricField(), Interp(new G4CMPTriLinearInterp(tli)) {;}

G4CMPMeshElectricField& 
G4CMPMeshElectricField::operator=(const G4CMPMeshElectricField &p) {
  if (this != &p) {				// Only copy if not self
    G4ElectricField::operator=(p);		// Call through to base
    *Interp = *p.Interp;
  }

  return *this;
}

G4CMPMeshElectricField::~G4CMPMeshElectricField() {
  delete Interp;
}

void G4CMPMeshElectricField::
BuildInterp(const std::vector<std::array<G4double,3> >& xyz,
	    const std::vector<G4double>& v,
	    const std::vector<std::array<G4int,4> >& tetra) {
  Interp->UseMesh(xyz, v, tetra);
}

void G4CMPMeshElectricField::BuildInterp(const G4String& EPotFileName,
					 G4double VScale) {
  if (G4CMPConfigManager::GetVerboseLevel() > 0) {
    G4cout << "G4CMPMeshElectricField::Constructor: Creating Electric Field " 
	   << EPotFileName;

    if (VScale != 1.) G4cout << " rescaled by " << VScale;
    G4cout << G4endl;
  }

  vector<array<G4double,4> > tempX;

  array<G4double,4> temp = {{ 0, 0, 0, 0 }};
  G4double x,y,z,v;

  G4double vmin=99999., vmax=-99999.;

  std::ifstream epotFile(EPotFileName);
  if (!epotFile.good()) {
    G4ExceptionDescription msg;
    msg << "Unable to open " << EPotFileName;
    G4Exception("G4CMPMeshElectricField::BuildInterp", "G4CMPEM001",
		FatalException, msg);
    return;
  }

  while (epotFile.good() && !epotFile.eof())
  {
    epotFile >> x >> y >> z >> v;
    temp[0] = x*m;
    temp[1] = y*m;
    temp[2] = z*m;
    temp[3] = v*volt * VScale;
    tempX.push_back(temp);

    if (temp[3]<vmin) vmin = temp[3];
    if (temp[3]>vmax) vmax = temp[3];
  }
  epotFile.close();

  if (G4CMPConfigManager::GetVerboseLevel() > 1) {
    G4cout << " Voltage from " << vmin/volt << " to " << vmax/volt << " V"
	   << G4endl;
  }

  std::sort(tempX.begin(),tempX.end(), vector_comp);

  vector<array<G4double, 3> > X(tempX.size(), {{0,0,0}});
  vector<G4double> V(tempX.size(),0);

  for (size_t ii = 0; ii < tempX.size(); ++ii)
  {
    X[ii][0] = tempX[ii][0];
    X[ii][1] = tempX[ii][1];
    X[ii][2] = tempX[ii][2];
    V[ii] = tempX[ii][3];
  }

  Interp->UseMesh(X, V);
}


void G4CMPMeshElectricField::GetFieldValue(const G4double Point[3],
				     G4double *Efield) const {
  G4ThreeVector InterpField = Interp->GetGrad(Point,true);	// No messages
  for (size_t i = 0; i < 3; ++i) {
    Efield[i] = 0.0;
    Efield[3+i] = -1 * InterpField[i];
  }
}


G4double G4CMPMeshElectricField::GetPotential(const G4double Point[3]) const {
  return Interp->GetValue(Point);		// Allow "outside hull" messages
}


G4bool G4CMPMeshElectricField::vector_comp(const array<G4double,4>& p1,
                                           const array<G4double,4>& p2) {
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
