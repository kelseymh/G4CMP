//  G4CMPInterpolator.cc
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20160610  Extracted from old G4CMPNR.cc
//  20160628  Active code moved from .hh to .cc

#include "G4CMPInterpolator.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>	/* To get abs(int) */
#include <limits>
#include <vector>
using namespace std;

/* The code in this section comes from Numerical Recipes III (Press et. al.)
   and is at some points modified slightly to suit our purposes */


// """""""""""""""" interp_1d.h METHODS (Numerical Respies III) """"""""""""""""
// constructor
G4CMPVInterpolator::G4CMPVInterpolator(const vector<double> &x,
				       const vector<double> &y, int m)
  : n(x.size()), mm(m), jsav(0), cor(0), xx(x), yy(y) {
  dj = min(1, (int)sqrt(sqrt((double)n)));
}  

// bin location
int G4CMPVInterpolator::locate(double x) {
  int ju,jm,jl;
  if (n < 2 || mm < 2 || mm > n) throw("locate size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  jl=0;
  ju=n-1;
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if ((x >= xx[jm]) == ascnd) jl=jm;
    else ju=jm;
  }
  cor = abs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return max(0,min(n-mm,jl-((mm-2)>>1)));
}

int G4CMPVInterpolator::hunt(double x)
{
  int jl=jsav, jm, ju, inc=1;
  if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  if (jl < 0 || jl > n-1) {
    jl=0;
    ju=n-1;
  } else {
    if ((x >= xx[jl]) == ascnd) {
      for (;;) {
	ju = jl + inc;
	if (ju >= n-1) { ju = n-1; break;}
	else if ((x < xx[ju]) == ascnd) break;
	else {
	  jl = ju;
	  inc += inc;
	}
      }
    } else {
      ju = jl;
      for (;;) {
	jl = jl - inc;
	if (jl <= 0) { jl = 0; break;}
	else if ((x >= xx[jl]) == ascnd) break;
	else {
	  ju = jl;
	  inc += inc;
	}
      }
    }
  }
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if ((x >= xx[jm]) == ascnd) jl=jm;
    else ju=jm;
  }
  cor = abs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return max(0,min(n-mm,jl-((mm-2)>>1)));
}
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

// >>>>>>>>>>>>>>>>> interp_linear.h from Numerical Recipes >>>>>>>>>>>>>>>>>>>
double G4CMPLinearInterp::rawinterp(int j, double x) {
  if (xx[j]==xx[j+1]) return yy[j];
  else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; interp_2d.h ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
G4CMPGridInterp& 
G4CMPGridInterp::operator=(const G4CMPGridInterp& oldBI) {
  if (this != &oldBI) {
    m = oldBI.m;
    n = oldBI.n;
    y = oldBI.y;
    x1terp = oldBI.x1terp;
    x2terp = oldBI.x2terp;
  }
  return *this;
}
    
double G4CMPGridInterp::interp(double x1p, double x2p) {
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
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
