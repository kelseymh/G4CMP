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
// 20190513  Provide support for 2D (e.g., axisymmetric) and 3D meshes.
// 20190919  BUG FIX:  2D project functions need 'break' in switch statements.
// 20200519  Move local "static" buffers to class for thread safety.
// 20210323  For 2D radial fields, need to manually protect rho < 0.

#include "G4CMPMeshElectricField.hh"
#include "G4CMPBiLinearInterp.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPTriLinearInterp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

using std::array;
using std::vector;


// Constructors

G4CMPMeshElectricField::
G4CMPMeshElectricField(const G4String& EPotFileName, G4double Vscale)
  : G4ElectricField(), Interp(0), xCoord(kUndefined), yCoord(kUndefined) {
  BuildInterp(EPotFileName, Vscale);
}

// Constructor for predefined 3D mesh table (or 2D projective mesh)

G4CMPMeshElectricField::
G4CMPMeshElectricField(const vector<array<G4double,3> >& xyz,
		       const vector<G4double>& v,
		       const vector<array<G4int,4> >& tetra,
		       EAxis xdim, EAxis ydim)
  : G4ElectricField(), Interp(0), xCoord(xdim), yCoord(ydim) {
  if (xdim == kUndefined) {	// Assume true 3D Cartesian mesh
    BuildInterp(xyz, v, tetra);
  } else {			// Projected 2D mesh with specified coordinates
    Interp = new G4CMPBiLinearInterp(xyz, v, tetra);
  }
}

// Constructor for predefined 2D mesh table with specified coordinates

G4CMPMeshElectricField::
G4CMPMeshElectricField(const std::vector<std::array<G4double,2> >& xy,
		       const std::vector<G4double>& v,
		       const std::vector<std::array<G4int,3> >& tetra,
		       EAxis xdim, EAxis ydim)
  : G4ElectricField(), Interp(0), xCoord(xdim), yCoord(ydim) {
  BuildInterp(xy, v, tetra);
}

// Copy previously constructed mesh interpolator

G4CMPMeshElectricField::G4CMPMeshElectricField(const G4CMPTriLinearInterp& tli)
  : G4ElectricField(), Interp(tli.Clone()), xCoord(kUndefined),
    yCoord(kUndefined) {;}

G4CMPMeshElectricField::G4CMPMeshElectricField(const G4CMPBiLinearInterp& bli,
					       EAxis xdim, EAxis ydim)
  : G4ElectricField(), Interp(bli.Clone()), xCoord(xdim), yCoord(ydim) {;}

G4CMPMeshElectricField::
G4CMPMeshElectricField(const G4CMPVMeshInterpolator* mesh,
		       EAxis xdim, EAxis ydim)
  : G4ElectricField(), Interp(mesh->Clone()), xCoord(xdim), yCoord(ydim) {;}

// Copy constructor and assignment operator

G4CMPMeshElectricField::G4CMPMeshElectricField(const G4CMPMeshElectricField &p)
  : G4ElectricField(p), Interp(p.Interp->Clone()), xCoord(p.xCoord),
    yCoord(p.yCoord) {;}

G4CMPMeshElectricField& 
G4CMPMeshElectricField::operator=(const G4CMPMeshElectricField &p) {
  if (this != &p) {				// Only copy if not self
    G4ElectricField::operator=(p);		// Call through to base

    if (Interp) delete Interp;
    Interp = p.Interp->Clone();
    xCoord = p.xCoord;
    yCoord = p.yCoord;
  }

  return *this;
}

G4CMPMeshElectricField::~G4CMPMeshElectricField() {
  delete Interp;
}


// Create 3D or 2D mesh interpolators from preconstructed mesh tables

void G4CMPMeshElectricField::
BuildInterp(const std::vector<std::array<G4double,3> >& xyz,
	    const std::vector<G4double>& v,
	    const std::vector<std::array<G4int,4> >& tetra) {
  if (!Interp) Interp = new G4CMPTriLinearInterp(xyz, v, tetra);
  else ((G4CMPTriLinearInterp*)Interp)->UseMesh(xyz, v, tetra);
}

void G4CMPMeshElectricField::
BuildInterp(const std::vector<std::array<G4double,2> >& xy,
	    const std::vector<G4double>& v,
	    const std::vector<std::array<G4int,3> >& tetra) {
  if (!Interp) Interp = new G4CMPBiLinearInterp(xy, v, tetra);
  else ((G4CMPBiLinearInterp*)Interp)->UseMesh(xy, v, tetra);
}


// Construct mesh from 3D input file

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

  while (epotFile.good() && !epotFile.eof()) {
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
 
  vector<array<G4double,3> > X(tempX.size(), {{0,0,0}});
  vector<G4double> V(tempX.size(),0);
  for (size_t ii = 0; ii < tempX.size(); ++ii)
  {
    X[ii][0] = tempX[ii][0];
    X[ii][1] = tempX[ii][1];
    X[ii][2] = tempX[ii][2];
    V[ii] = tempX[ii][3];
  }
 
  if (Interp) delete Interp;
  Interp = new G4CMPTriLinearInterp(X, V);
}


// Use mesh interpolater to evaluate electric field

void G4CMPMeshElectricField::GetFieldValue(const G4double Point[3],
					   G4double *BEfield) const {
  G4ThreeVector InterpField;
  if (xCoord == kUndefined) {		// Three dimensions
    InterpField = Interp->GetGrad(Point,true);
  } else {				// Two dimensions
    G4double proj[2] = { 0.,0. };
    Project2D(Point, proj);
    InterpField = Interp->GetGrad(proj, true);
    Expand2Dat(Point, InterpField);
  }

  for (size_t i = 0; i < 3; ++i) {
    BEfield[i] = 0.0;
    BEfield[3+i] = -1 * InterpField[i];
  }
}

G4double G4CMPMeshElectricField::GetPotential(const G4double Point[3]) const {
  if (xCoord == kUndefined) {		// Three dimensions
    return Interp->GetValue(Point);
  } else {				// Two dimensions
    G4double proj[2] = { 0.,0. };
    Project2D(Point, proj);
    return Interp->GetValue(proj);
  }
}


// Convert between 3D and 2D coordinates for projected meshes

namespace {
  const char* AxisName(EAxis axis) {	// Convenience function for debugging
    return (axis==kXAxis    ? "kXAxis"    : axis==kYAxis ? "kYAxis" :
	    axis==kZAxis    ? "kZAxis"    : axis==kRho   ? "kRho" :
	    axis==kRadial3D ? "kRadial3D" : axis==kPhi   ? "kPhi" :
	    "UNKNOWN");
  }
}

void G4CMPMeshElectricField::Project2D(const G4double Point[3],
				       G4double Project[2]) const {
  pos_.set(Point[0],Point[1],Point[2]);

  Project[0] = Project[1] = 0.;
  
  switch (xCoord) {
  case kXAxis:    Project[0] = Point[0];  break;
  case kYAxis:    Project[0] = Point[1];  break;
  case kZAxis:    Project[0] = Point[2];  break;
  case kRho:      Project[0] = pos_.rho(); break;
  case kRadial3D: Project[0] = pos_.r();   break;
  case kPhi:	  Project[0] = pos_.phi(); break;
  default: ;
  }
  
  switch (yCoord) {
  case kXAxis:    Project[1] = Point[0];  break;
  case kYAxis:    Project[1] = Point[1];  break;
  case kZAxis:    Project[1] = Point[2];  break;
  case kRho:      Project[1] = pos_.rho(); break;
  case kRadial3D: Project[1] = pos_.r();   break;
  case kPhi:	  Project[1] = pos_.phi(); break;
  default: ;
  }

  if (G4CMPConfigManager::GetVerboseLevel() > 2) {
    G4cout << "Project2D Point " << pos_ << " onto axes " << AxisName(xCoord)
	   << " " << AxisName(yCoord) << " : (" << Project[0] << ","
	   << Project[1] << ")" << G4endl;
  }
}

void G4CMPMeshElectricField::Expand2Dat(const G4double Point[3],
					G4ThreeVector& Efield) const {
  pos_.set(Point[0],Point[1],Point[2]);

  G4double xval = Efield.x(), yval = Efield.y();

  // Cartesian projections are easy and direct
  if (xCoord == kXAxis && yCoord == kYAxis) Efield.set(xval, yval, 0.);
  if (xCoord == kXAxis && yCoord == kZAxis) Efield.set(xval, 0., yval);
  if (xCoord == kYAxis && yCoord == kZAxis) Efield.set(0., xval, yval);

  // Radial field (e.g., spherical electrode) is easy
  // NOTE: CLHEP doesn't like R<0; flip angles manually to compensate
  if (xCoord == kRadial3D) {
    Efield.setRThetaPhi(fabs(xval), xval<0?pi-pos_.theta():pos_.theta(),
			pos_.phi()+(xval<0?pi:0.));
  }

  // Cylindrical field around Z axis is easy
  // NOTE: CLHEP doesn't like R<0; flip angle manually to compensate
  if (xCoord == kRho && yCoord == kZAxis) {
    Efield.setRhoPhiZ(fabs(xval), pos_.phi()+(xval<0?pi:0.), yval);
  }

  // Cylindrical fields around other axes are more complicated
  if (xCoord == kRho && yCoord == kXAxis) {
    G4double phiYZ = atan2(pos_.z(), pos_.y());
    Efield.set(yval, xval*cos(phiYZ), xval*sin(phiYZ));
  }
  
  if (xCoord == kRho && yCoord == kYAxis) {
    G4double phiZX = atan2(pos_.x(), pos_.z());
    Efield.set(xval*sin(phiZX), yval, xval*cos(phiZX));
  }
}


// Sorting function for reading meash points from 3D input file

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
