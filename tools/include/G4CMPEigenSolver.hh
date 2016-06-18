//  G4CMPEigenSolver.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef _G4CMPEigenSolver_hh
#define _G4CMPEigenSolver_hh

/* The code in this section comes from Numerical Recipes III (Press et. al.)
   and is at some points modified slightly to suit our purposes */

// some system #include's we'll need
#include "matrix.hh"
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

  G4CMPEigenSolver(const matrix<double> &a, bool yesvec=true) { init(a, yesvec); }

  void init(const matrix<double> &a, bool yesvec=true) {
    n = a.nrows();
    z = a;
    d.resize(n,0.);
    e.resize(n,0.);
    yesvecs = yesvec;

    tred2();
    tqli();
    sort();
  }

  G4CMPEigenSolver(const vector<double> &dd, const vector<double> &ee, bool yesvec=true) {
    init(dd, ee, yesvec);
  }

  void init(const vector<double> &dd, const vector<double> &ee, bool yesvec=true) {
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

#endif	/* G4CMPEigenSolver_hh */
