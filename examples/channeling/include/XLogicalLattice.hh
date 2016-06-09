/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/XLogicalLattice.hh
/// \brief Definition of the XLogicalLattice class
//
// $Id: 54c8177a418708cbc01443ad389dade2f2cedfd1 $
//
#ifndef XLogicalLattice_h
#define XLogicalLattice_h

#define MAXRES 322  //maximum one dimensional map resolution

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"


using namespace std;

class XLogicalLattice{

private:
  
  int fVresTheta; //velocity  map theta resolution (inclination)
  int fVresPhi;   //velocity  map phi resolution  (azimuth)
  int fDresTheta; //direction map theta resn
  int fDresPhi;   //direction map phi resn 

  double fMap[3][MAXRES][MAXRES];  //map for group velocity scalars
  double fN_map[3][MAXRES][MAXRES][3];//normalized map containing group V direction unit vectors

  double fA;       //Scaling constant for Anh.Dec. mean free path
  double fB;       //Scaling constant for Iso.Scat. mean free path
  double fDosL;    //Density of states for L-phonons
  double fDosST;   //Density of states for ST-phonons
  double fDosFT;    //Density of states for FT-phonons
  double fBeta, fGamma, fLambda, fMu; //dynamical constants for material


  ifstream fMapFile;
  int fThetaRes, fPhiRes;

public:


  
  XLogicalLattice();
  ~XLogicalLattice();

  void SetDynamicalConstants(double, double, double, double);
  void SetScatteringConstant(G4double);
  void SetAnhDecConstant(G4double);
  void SetLDOS(double);
  void SetSTDOS(double);
  void SetFTDOS(double);

    
    
  double GetBeta();
  double GetGamma();
  double GetLambda();
  double GetMu();
  G4double GetScatteringConstant();
  G4double GetAnhDecConstant();
  double GetLDOS();
  double GetSTDOS();
  double GetFTDOS();

    
  bool LoadMap(int, int, int, string);
  bool Load_NMap(int, int, int, string);
  double MapKtoV(int, G4ThreeVector);   //Get full group velocity vector
  G4ThreeVector MapKtoVDir(int, G4ThreeVector);//Get normalized group velocity direction so that normalisatioon does not have to be done at run time

};

#endif
