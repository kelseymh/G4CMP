//  G4CMPInterpolator.h
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160610  Extracted from old G4CMPNR.hh

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
  const double *xx, *yy;

  G4CMPVInterpolator(const vector<double> &x, const double *y, int m);
    
  double interp(double x) {
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
    : G4CMPVInterpolator(xv,&yv[0],2)  {}

  double rawinterp(int j, double x) {
    if (xx[j]==xx[j+1]) return yy[j];
    else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
  }
};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; interp_2d.h ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
struct G4CMPBiLinearInterp {
  int m,n;
  matrix<double> y;
  G4CMPLinearInterp x1terp, x2terp;
  
  G4CMPBiLinearInterp(const vector<double> &x1v, const vector<double> &x2v,
		      const matrix<double> &ym)
    : m(x1v.size()), n(x2v.size()), y(ym), x1terp(x1v,x1v), x2terp(x2v,x2v) {}
    
  // adding this method seems to be required to make things run...
  // ... I am not sure why or why having it empty works, but it seems to
  G4CMPBiLinearInterp& operator=(const G4CMPBiLinearInterp& oldBI) {
    // Proper content added by M. Kelsey
    m = oldBI.m;
    n = oldBI.n;
    y = oldBI.y;
    x1terp = oldBI.x1terp;
    x2terp = oldBI.x2terp;
    return *this;
  }
    
  double interp(double x1p, double x2p) {
    int i,j;
    double yy, t, u;
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
    t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
    u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
    yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j]
      + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
    return yy;
  }
};
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#endif	/* G4CMPInterpolator_hh */
