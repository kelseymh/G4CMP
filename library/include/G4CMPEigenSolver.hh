//  G4CMPEigenSolver.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef _G4CMPEigenSolver_hh
#define _G4CMPEigenSolver_hh

/* The code in this section comes from Numerical Recipes III (Press et. al.)
   and is at some points modified slightly to suit our purposes */

// some system #include's we'll need
#include "G4CMPMatrix.hh"
#include <vector>
using G4CMP::matrix;
using std::vector;

// ................. eigen_sym.h from Numerical Recipes III ...................
struct G4CMPEigenSolver {
  size_t n;
  matrix<double> z;
  vector<double> d,e;
  bool yesvecs;

  G4CMPEigenSolver() : n(0), yesvecs(false) {;}

  G4CMPEigenSolver(const matrix<double> &a, bool yesvec=true)
    : n(a.rows()), z(a), d(a.rows(),0.), e(a.rows(),0.), yesvecs(yesvec) {
    tred2();
    tqli();
    sort();
  }

  G4CMPEigenSolver(const vector<double> &dd, const vector<double> &ee,
		   bool yesvec=true)
    : n(dd.size()), z(dd.size(),dd.size(),0.), d(dd), e(ee), yesvecs(yesvec) {
    for (size_t i=0;i<n;i++) z[i][i]=1.0;
    tqli();
    sort();
  }

  // Reusable with matrix constructor above
  void setup(const matrix<double> &a, bool yesvec=true) {
    n = a.rows();
    z = a;
    d.clear(); d.resize(n,0.);
    e.clear(); e.resize(n,0.);
    yesvecs = yesvec;
    tred2();
    tqli();
    sort();
  }

  void sort();
  void tred2();
  void tqli();
  double pythag(double a, double b);
};
// ............................................................................

#endif	/* G4CMPEigenSolver_hh */
