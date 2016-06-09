//  G4CMPEigenSolver.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef _G4CMPEigenSolver_h
#define _G4CMPEigenSolver_h

/* The code in this section comes from Numerical Recipes III (Press et. al.)
   and is at some points modified slightly to suit our purposes */

// some system #include's we'll need
#include "matrix.hh"
#include <vector>

using namespace std;

/*****
// exception handling
#ifndef _USENRERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
	char *message;
	char *file;
	int line;
	NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
           err.message, err.file, err.line);
	exit(1);
}
#endif
*****/

// @@@@@@@@@@@@@@@@@@@@@@@@@@@  Matrix Class @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

typedef vector<double> VecDoub;
typedef matrix<int> MatInt;
typedef matrix<double> MatDoub;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


// ................. eigen_sym.h from Numerical Recipes III ...................
struct G4CMPEigenSolver {
  size_t n;
  MatDoub z;
  VecDoub d,e;
  bool yesvecs;

  G4CMPEigenSolver() : n(0), yesvecs(false) {;}

  G4CMPEigenSolver(const MatDoub &a, bool yesvec=true) { init(a, yesvec); }

  void init(const MatDoub &a, bool yesvec=true) {
    n = a.nrows();
    z = a;
    d.resize(n,0.);
    e.resize(n,0.);
    yesvecs = yesvec;

    tred2();
    tqli();
    sort();
  }

  G4CMPEigenSolver(const VecDoub &dd, const VecDoub &ee, bool yesvec=true) {
    init(dd, ee, yesvec);
  }

  void init(const VecDoub &dd, const VecDoub &ee, bool yesvec=true) {
    n = dd.size();
    d = dd;
    e = ee;
    yesvecs = yesvec;

    z.resize(n,n,0.0);
    for (int i=0;i<n;i++) z[i][i]=1.0;
    tqli();
    sort();
  }

  void sort();
  void tred2();
  void tqli();
  double pythag(double a, double b);
};
// ............................................................................


// """""""""""""""""" interp_1d.h from Numerical Respies III """"""""""""""""""
struct G4CMPVInterpolator {
  int n, mm, jsav, cor, dj;
  const double *xx, *yy;

  G4CMPVInterpolator(const VecDoub &x, const double *y, int m);
    
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
  G4CMPLinearInterp(const VecDoub &xv, const VecDoub &yv)
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
  MatDoub y;
  G4CMPLinearInterp x1terp, x2terp;
  
  G4CMPBiLinearInterp(const VecDoub &x1v, const VecDoub &x2v, const MatDoub &ym)
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

#endif	/* G4CMPEigenSolver.h */
