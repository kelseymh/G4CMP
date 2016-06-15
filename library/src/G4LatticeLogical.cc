/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file materials/src/G4LatticeLogical.cc
/// \brief Implementation of the G4LatticeLogical class
//
// $Id$
//
// 20140218  Add new charge-carrier parameters to output
// 20140306  Allow valley filling using Euler angles directly
// 20140318  Compute electron mass scalar (Herring-Vogt) from tensor
// 20140324  Include inverse mass-ratio tensor
// 20140408  Add valley momentum calculations
// 20140425  Add "effective mass" calculation for electrons
// 20150601  Add mapping from electron velocity back to momentum
// 20160517  Add basis vectors for lattice, to use with Miller orientation
// 20160520  Add reporting function to format valley Euler angles
// 20160614  Add elasticity tensors and density (set from G4Material) 

#include "G4LatticeLogical.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeLogical::G4LatticeLogical(const G4String& name)
  : verboseLevel(0), fName(name),
    fDensity(0.), fElasticity{}, fElReduced{}, fHasElasticity(false),
    fVresTheta(0), fVresPhi(0), fDresTheta(0), fDresPhi(0),
    fA(0), fB(0), fLDOS(0), fSTDOS(0), fFTDOS(0),
    fBeta(0), fGamma(0), fLambda(0), fMu(0),
    fVSound(0.), fL0_e(0.), fL0_h(0.), 
    mElectron(electron_mass_c2/c_squared),
    fHoleMass(mElectron), fElectronMass(mElectron),
    fMassTensor(G4Rep3x3(mElectron,0.,0.,0.,mElectron,0.,0.,0.,mElectron)),
    fMassInverse(G4Rep3x3(1/mElectron,0.,0.,0.,1/mElectron,0.,0.,0.,1/mElectron)),
    fIVField(0.), fIVRate(0.), fIVExponent(0.) {
  for (G4int i=0; i<3; i++) {
    for (G4int j=0; j<MAXRES; j++) {
      for (G4int k=0; k<MAXRES; k++) {
	fMap[i][j][k] = 0.;
	fN_map[i][j][k].set(0.,0.,0.);
      }
    }
  }
  SetElasticityCubic(0.,0.,0.);		// Fill elasticity tensors with zeros;
  fHasElasticity = false;
}

G4LatticeLogical::~G4LatticeLogical() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
//Complete basis vectors: right-handed, possibly orthonormal
/////////////////////////////////////////////////////////////
void G4LatticeLogical::SetBasis() {
  static const G4ThreeVector origin(0.,0.,0.);
  if (fBasis[0].isNear(origin,1e-6)) fBasis[0].set(1.,0.,0.);
  if (fBasis[1].isNear(origin,1e-6)) fBasis[1].set(0.,1.,0.);
  if (fBasis[2].isNear(origin,1e-6)) fBasis[2] = fBasis[0].cross(fBasis[1]);

  // Ensure that all basis vectors are unit
  fBasis[0].setMag(1.);
  fBasis[1].setMag(1.);
  fBasis[2].setMag(1.);

  if (fBasis[0].cross(fBasis[1]).dot(fBasis[2]) < 0.) {
    G4cerr << "ERROR G4LatticeLogical has a left-handed basis!" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/////////////////////////////////////////////////////////////
// Elasticity matrix, with cubic symmetries
/////////////////////////////////////////////////////////////
void 
G4LatticeLogical::SetElasticityCubic(G4double C11, G4double C12, G4double C44) {
  if (verboseLevel) {
    G4cout << "G4LatticeLogical[" << fName << "]::SetElasticityCubic "
	   << C11 << " " << C12 << " " << C44 << G4endl;
  }

  // Reduced elasticity tensor is block-symmetric 6x6 array
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      fElReduced[i][j] = (i==j) ? C11 : C12;
    }
  }

  for (int i=3; i<6; i++) {
    fElReduced[i][i] = C44;
  }

  // Unpack reduced elasticity tensor into full four-dimensional Cijkl
  G4int rn1[6][2] = { };
  rn1[1][0] = rn1[1][1] = rn1[3][0] = rn1[5][1] = 1;
  rn1[2][0] = rn1[2][1] = rn1[3][1] = rn1[4][0] = 2;

  G4int rn2[2][2] = { };
  rn2[0][1] = rn2[1][0] = 1;

  for(int l=0; l<2; l++) {
    for(int k=0; k<2; k++) {
      for(int j=0; j<6; j++) {
	for(int i=0; i<6; i++) {
	  fElasticity[rn1[i][rn2[k][0]]][rn1[i][rn2[k][1]]]
	             [rn1[j][rn2[l][0]]][rn1[j][rn2[l][1]]] = fElReduced[i][j];
	}
      }
    }
  }

  fHasElasticity = true;	// Flag use of tensors is safe
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

///////////////////////////////////////////
//Load map of group velocity scalars (m/s)
////////////////////////////////////////////
G4bool G4LatticeLogical::LoadMap(G4int tRes, G4int pRes,
				 G4int polarizationState, G4String map) {
  if (tRes>MAXRES || pRes>MAXRES) {
    G4cerr << "G4LatticeLogical::LoadMap exceeds maximum resolution of "
           << MAXRES << " by " << MAXRES << ". terminating." << G4endl;
    return false; 		//terminate if resolution out of bounds.
  }

  std::ifstream fMapFile(map.data());
  if (!fMapFile.is_open()) return false;

  G4double vgrp = 0.;
  for (G4int theta = 0; theta<tRes; theta++) {
    for (G4int phi = 0; phi<pRes; phi++) {
      fMapFile >> vgrp;
      fMap[polarizationState][theta][phi] = vgrp * (m/s);
    }
  }

  if (verboseLevel) {
    G4cout << "\nG4LatticeLogical::LoadMap(" << map << ") successful"
	   << " (Vg scalars " << tRes << " x " << pRes << " for polarization "
	   << polarizationState << ")." << G4endl;
  }

  fVresTheta=tRes; //store map dimensions
  fVresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


////////////////////////////////////
//Load map of group velocity unit vectors
///////////////////////////////////
G4bool G4LatticeLogical::Load_NMap(G4int tRes, G4int pRes,
				   G4int polarizationState, G4String map) {
  if (tRes>MAXRES || pRes>MAXRES) {
    G4cerr << "G4LatticeLogical::LoadMap exceeds maximum resolution of "
           << MAXRES << " by " << MAXRES << ". terminating." << G4endl;
    return false; 		//terminate if resolution out of bounds.
  }

  std::ifstream fMapFile(map.data());
  if(!fMapFile.is_open()) return false;

  G4double x,y,z;	// Buffers to read coordinates from file
  G4ThreeVector dir;
  for (G4int theta = 0; theta<tRes; theta++) {
    for (G4int phi = 0; phi<pRes; phi++) {
      fMapFile >> x >> y >> z;
      dir.set(x,y,z);
      fN_map[polarizationState][theta][phi] = dir.unit();	// Enforce unity
    }
  }

  if (verboseLevel) {
    G4cout << "\nG4LatticeLogical::Load_NMap(" << map << ") successful"
	   << " (Vdir " << tRes << " x " << pRes << " for polarization "
	   << polarizationState << ")." << G4endl;
  }

  fDresTheta=tRes; //store map dimensions
  fDresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon velocity in m/s

G4double G4LatticeLogical::MapKtoV(G4int polarizationState,
				   const G4ThreeVector& k) const {
  G4double theta, phi, tRes, pRes;

  tRes=pi/fVresTheta;
  pRes=twopi/fVresPhi;
  
  theta=k.getTheta();
  phi=k.getPhi();

  if(phi<0) phi = phi + twopi;
  if(theta>pi) theta=theta-pi;

  G4double Vg = fMap[polarizationState][int(theta/tRes)][int(phi/pRes)];

  if (Vg == 0) {
    G4cerr << "Found v=0 for polarization "<< polarizationState
	   << " theta " << theta << " phi " << phi
	   << " translating to map coords"
	   << " theta " << int(theta/tRes) << " phi " << int(phi/pRes)
	   << G4endl;
  }

  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::MapKtoV theta,phi=" << theta << " " << phi
	   << " : ith,iph " << int(theta/tRes) << " " << int(phi/pRes)
	   << " : V " << Vg << G4endl;
  }

  return Vg;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon propagation direction as dimensionless unit vector

G4ThreeVector G4LatticeLogical::MapKtoVDir(G4int polarizationState,
					   const G4ThreeVector& k) const {  
  G4double theta, phi, tRes, pRes;

  tRes=pi/(fDresTheta-1);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fDresPhi-1);

  theta=k.getTheta();
  phi=k.getPhi(); 

  if(theta>pi) theta=theta-pi;
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;
  if(phi<0) phi = phi + 2*pi;

  G4int iTheta = int(theta/tRes+0.5);
  G4int iPhi = int(phi/pRes+0.5);

  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::MapKtoVDir theta,phi=" << theta << " " << phi
	   << " : ith,iph " << iTheta << " " << iPhi
	   << " : dir " << fN_map[polarizationState][iTheta][iPhi] << G4endl;
  }

  return fN_map[polarizationState][iTheta][iPhi];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Convert electron momentum to valley velocity, wavevector, and HV vector

G4ThreeVector 
G4LatticeLogical::MapPtoV_el(G4int ivalley, const G4ThreeVector& p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoV_el " << ivalley << " " << p_e
	   << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return vToN.inverse()*(GetMInvTensor()*(vToN*p_e/c_light));
}

G4ThreeVector 
G4LatticeLogical::MapV_elToP(G4int ivalley, const G4ThreeVector& v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToP " << ivalley << " " << v_e
	   << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return vToN.inverse()*(GetMassTensor()*(vToN*v_e*c_light));
}

G4ThreeVector
G4LatticeLogical::MapV_elToK_HV(G4int ivalley, const G4ThreeVector &v_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToK_HV " << ivalley << " " << v_e
     << G4endl;

  const G4RotationMatrix& vToN = GetValley(ivalley);
  return GetSqrtInvTensor()*GetMassTensor()*vToN*v_e/hbar_Planck;
}

G4ThreeVector 
G4LatticeLogical::MapPtoK_valley(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoK " << ivalley << " " << p_e
	   << G4endl;

  p_e /= hbarc;					// Convert to wavevector
  return p_e.transform(GetValley(ivalley));	// Rotate into valley frame
}

G4ThreeVector 
G4LatticeLogical::MapPtoK_HV(G4int ivalley, G4ThreeVector p_e) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoK_HV " << ivalley << " " << p_e
	   << G4endl;

  p_e.transform(GetValley(ivalley));		// Rotate into valley frame
  return GetSqrtInvTensor() * p_e/hbarc;	// Herring-Vogt transformation
}

G4ThreeVector 
G4LatticeLogical::MapK_HVtoK_valley(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoK_valley " << ivalley << " " << k_HV
	   << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley
  return k_HV;
}

G4ThreeVector
G4LatticeLogical::MapK_HVtoK(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoK " << ivalley << " " << k_HV
     << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley
  k_HV.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  return k_HV;
}

G4ThreeVector
G4LatticeLogical::MapK_HVtoV_el(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoV_el " << ivalley << " " << k_HV
     << G4endl;

  k_HV *= hbar_Planck; // k_HV to p_HV
  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley
  k_HV *= GetMInvTensor();			// From p_valley to v_valley
  k_HV.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  return k_HV;
}

G4ThreeVector 
G4LatticeLogical::MapK_HVtoP(G4int ivalley, G4ThreeVector k_HV) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_HVtoP " << ivalley << " " << k_HV
	   << G4endl;

  k_HV *= GetSqrtTensor();			// From Herring-Vogt to valley 
  k_HV.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  k_HV *= hbarc;			// Convert wavevector to momentum
  return k_HV;
}

G4ThreeVector 
G4LatticeLogical::MapK_valleyToP(G4int ivalley, G4ThreeVector k) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapK_valleyToP " << ivalley << " " << k
	   << G4endl;

  k.transform(GetValley(ivalley).inverse());	// Rotate out of valley frame
  k *= hbarc;				// Convert wavevector to momentum
  return k;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Apply energy-momentum relationship for electron transport

G4double  
G4LatticeLogical::MapPtoEkin(G4int iv, G4ThreeVector p) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapPtoEkin " << iv << " " << p << G4endl;

  p.transform(GetValley(iv));			// Rotate to valley frame

  // Compute kinetic energy component by component, then sum
  return (0.5/c_squared) * (p.x()*p.x()*fMassInverse.xx() +
			    p.y()*p.y()*fMassInverse.yy() +
			    p.z()*p.z()*fMassInverse.zz());
}

G4double
G4LatticeLogical::MapV_elToEkin(G4int iv, G4ThreeVector v) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::MapV_elToEkin " << iv << " " << v << G4endl;

  v.transform(GetValley(iv));			// Rotate to valley frame

  // Compute kinetic energy component by component, then sum
  return 0.5 * (v.x()*v.x()*fMassTensor.xx() +
          v.y()*v.y()*fMassTensor.yy() +
          v.z()*v.z()*fMassTensor.zz());
}

// Compute effective "scalar" electron mass to match energy/momentum relation

G4double 
G4LatticeLogical::GetElectronEffectiveMass(G4int iv,
					   const G4ThreeVector& p) const {
  if (verboseLevel>1)
    G4cout << "G4LatticeLogical::GetElectronEffectiveMass " << iv
	   << " " << p << G4endl;

  return 0.5*p.mag2()/c_squared/MapPtoEkin(iv,p);	// Non-relativistic
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Store electron mass tensor using diagonal elements

void G4LatticeLogical::SetMassTensor(G4double mXX, G4double mYY, G4double mZZ) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::SetMassTensor " << mXX << " " << mYY
	   << " " << mZZ << " *m_e" << G4endl;
  }

  // NOTE:  Use of G4RotationMatrix not appropriate here, as matrix is
  //        not normalized.  But CLHEP/Matrix not available in GEANT4.
  fMassTensor.set(G4Rep3x3(mXX*mElectron, 0., 0.,
			   0., mYY*mElectron, 0.,
			   0., 0., mZZ*mElectron));

  FillMassInfo();
}

void G4LatticeLogical::SetMassTensor(const G4RotationMatrix& etens) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::SetMassTensor " << etens << G4endl;
  }

  // Check if mass tensor already has electron mass, or is just coefficients
  G4bool hasEmass = (etens.xx()/mElectron > 1e-3 ||
		     etens.yy()/mElectron > 1e-3 ||
		     etens.zz()/mElectron > 1e-3);
  G4double mscale = hasEmass ? 1. : mElectron;

  // NOTE:  Use of G4RotationMatrix not appropriate here, as matrix is
  //        not normalized.  But CLHEP/Matrix not available in GEANT4.
  fMassTensor.set(G4Rep3x3(etens.xx()*mscale, 0., 0.,
			   0., etens.yy()*mscale, 0.,
			   0., 0., etens.zz()*mscale));

  FillMassInfo();
}

// Compute derived quantities from user-input mass tensor

void G4LatticeLogical::FillMassInfo() {
  // Herring-Vogt scalar mass of electron, used for sound speed calcs  
  fElectronMass = 3. / ( 1./fMassTensor.xx() + 1./fMassTensor.yy()
			 + 1./fMassTensor.zz() );  

  // 1/m mass tensor used for k and v calculations in valley coordinates
  fMassInverse.set(G4Rep3x3(1./fMassTensor.xx(), 0., 0.,
			    0., 1./fMassTensor.yy(), 0.,
			    0., 0., 1./fMassTensor.zz()));

  // Mass ratio tensor used for scattering and field calculations
  fMassRatioSqrt.set(G4Rep3x3(sqrt(fMassTensor.xx()/fElectronMass), 0., 0.,
			      0., sqrt(fMassTensor.yy()/fElectronMass), 0.,
			      0., 0., sqrt(fMassTensor.zz()/fElectronMass)));

  fMInvRatioSqrt.set(G4Rep3x3(1./fMassRatioSqrt.xx(), 0., 0.,
			      0., 1./fMassRatioSqrt.yy(), 0.,
			      0., 0., 1./fMassRatioSqrt.zz()));
}

// Store drifting-electron valley using Euler angles

void G4LatticeLogical::AddValley(G4double phi, G4double theta, G4double psi) {
  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::AddValley " << psi << " " << theta
	   << " " << psi << " rad" << G4endl;
  }

  // Extend vector first, then fill last value, to reduce temporaries
  fValley.resize(fValley.size()+1);
  fValley.back().set(phi,theta,psi);
}

// Transform for drifting-electron valleys in momentum space

const G4RotationMatrix& G4LatticeLogical::GetValley(G4int iv) const {
  if (verboseLevel>1) G4cout << "G4LatticeLogical::GetValley " << iv << G4endl;

  if (iv >=0 && iv < (G4int)NumberOfValleys()) return fValley[iv];

  if (verboseLevel)
    G4cerr << "G4LatticeLogical ERROR: No such valley " << iv << G4endl;
  return G4RotationMatrix::IDENTITY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump structure in format compatible with reading back

void G4LatticeLogical::Dump(std::ostream& os) const {
  os << "# " << fName << " crystal lattice parameters" << std::endl;
  if (fHasElasticity) {		// TEMPORARY: Assume cubic only!
    os << "cubic " << fElReduced[1][1]
       << " " << fElReduced[1][2]
       << " " << fElReduced[4][4] << std::endl;
  }

  for (size_t i=0; i<3; i++) {
    os << "basis " << fBasis[0].x() << " " << fBasis[0].y()
       << " " << fBasis[0].z() << std::endl;
  }

  os << "# Phonon propagation parameters"
     << "\ndyn " << fBeta << " " << fGamma << " " << fLambda << " " << fMu
     << "\nscat " << fB << " decay " << fA
     << "\nLDOS " << fLDOS << " STDOS " << fSTDOS << " FTDOS " << fFTDOS
     << std::endl;

  os << "# Charge carrier propagation parameters"
     << "\nhmass " << fHoleMass/mElectron
     << "\nemass " << fMassTensor.xx()/mElectron
     << " " << fMassTensor.yy()/mElectron
     << " " << fMassTensor.zz()/mElectron << std::endl;

  os << "# Inverse mass tensor: " << fMassInverse.xx()*mElectron
     << " " << fMassInverse.yy()*mElectron
     << " " << fMassInverse.zz()*mElectron
     << " * 1/m(electron)" << std::endl
     << "# Herring-Vogt scalar mass: " << fElectronMass/mElectron << std::endl
     << "# sqrt(tensor/scalor): " << fMassRatioSqrt.xx()
     << " " << fMassRatioSqrt.yy()
     << " " << fMassRatioSqrt.zz()
     << std::endl;

  for (size_t i=0; i<NumberOfValleys(); i++) {
    DumpValley(os, i);
  }

  os << "# Intervalley scattering parameters"
     << "\nivField " << fIVField << "\t# V/m"
     << "\nivRate " << fIVRate/s << "\t# s"
     << "\nivPower" << fIVExponent << std::endl;

  os << "# Phonon wavevector/velocity maps" << std::endl;
  Dump_NMap(os, 0, "LVec.ssv");
  Dump_NMap(os, 1, "FTVec.ssv");
  Dump_NMap(os, 2, "STVec.ssv");

  DumpMap(os, 0, "L.ssv");
  DumpMap(os, 1, "FT.ssv");
  DumpMap(os, 2, "ST.ssv");
}

void G4LatticeLogical::DumpMap(std::ostream& os, G4int pol,
			       const G4String& name) const {
  os << "VG " << name << " " << (pol==0?"L":pol==1?"FT":pol==2?"ST":"??")
     << " " << fVresTheta << " " << fVresPhi << std::endl;

  if (verboseLevel) {
    for (G4int iTheta=0; iTheta<fVresTheta; iTheta++) {
      for (G4int iPhi=0; iPhi<fVresPhi; iPhi++) {
	os << fMap[pol][iTheta][iPhi] << std::endl;
      }
    }
  }
}

void G4LatticeLogical::Dump_NMap(std::ostream& os, G4int pol,
				 const G4String& name) const {
  os << "VDir " << name << " " << (pol==0?"L":pol==1?"FT":pol==2?"ST":"??")
     << " " << fDresTheta << " " << fDresPhi << std::endl;

  if (verboseLevel) {
    for (G4int iTheta=0; iTheta<fDresTheta; iTheta++) {
      for (G4int iPhi=0; iPhi<fDresPhi; iPhi++) {
	os << fN_map[pol][iTheta][iPhi].x()
	   << " " << fN_map[pol][iTheta][iPhi].y()
	   << " " << fN_map[pol][iTheta][iPhi].z()
	   << std::endl;
      }
    }
  }
}

// Print out Euler angles of requested valley

void G4LatticeLogical::DumpValley(std::ostream& os, G4int iv) const {
  if (iv < 0 || iv >= NumberOfValleys()) return;

  os << "valley " << fValley[iv].phi()/deg
     << " " << fValley[iv].theta()/deg
     << " " << fValley[iv].psi()/deg
     << " deg" << std::endl;
}
