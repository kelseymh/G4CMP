//  G4CMPInterpolator.h
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160610  Extracted from old G4CMPNR.hh
//  20160628  Active code moved from .hh to .cc
//  20181010  Address compiler warnings; define virtual dtors

#ifndef _G4CMPInterpolator_hh
#define _G4CMPInterpolator_hh

#include "G4CMPMatrix.hh"
#include <vector>
using G4CMP::matrix;
using std::vector;

/* The code in this section comes from Numerical Recipes III (Press et. al.)
   and is at some points modified slightly to suit our purposes */

// """""""""""""""""" interp_1d.h from Numerical Respies III """"""""""""""""""
struct G4CMPVInterpolator {
  int n, mm, jsav, cor, dj;
  vector<double> xx, yy;

  G4CMPVInterpolator(const vector<double>& x, const vector<double>& y, int m);
  virtual ~G4CMPVInterpolator() {;}

  virtual double interp(double x) {
    int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }
    
  int locate(double x);
  int hunt(double x);
  
  virtual double rawinterp(int jlo, double x) = 0;
};
// """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

// >>>>>>>>>>>>>>>>> interp_linear.h from Numerical Recipes >>>>>>>>>>>>>>>>>>>
struct G4CMPLinearInterp : G4CMPVInterpolator {
  G4CMPLinearInterp(const vector<double> &xv, const vector<double> &yv)
    : G4CMPVInterpolator(xv,yv,2)  {}
  virtual ~G4CMPLinearInterp() {;}

  virtual double rawinterp(int j, double x);
};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; interp_2d.h ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
struct G4CMPGridInterp {
  int m,n;
  matrix<double> y;
  G4CMPLinearInterp x1terp, x2terp;
  
  G4CMPGridInterp(const vector<double> &x1v, const vector<double> &x2v,
		      const matrix<double> &ym)
    : m(x1v.size()), n(x2v.size()), y(ym), x1terp(x1v,x1v), x2terp(x2v,x2v) {;}
  virtual ~G4CMPGridInterp() {;}

  G4CMPGridInterp& operator=(const G4CMPGridInterp& oldBI);
    
  virtual double interp(double x1p, double x2p);
};
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#endif	/* G4CMPInterpolator_hh */
