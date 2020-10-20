// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                        --- FanoBinomial ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// 20201019  Michael Kelsey (TAMU)
// =======================================================================

#include "G4CMPFanoBinomial.hh"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Random/RandBinomial.h"
#include <algorithm>	// for min() and max()
#include <cmath>	// for exp()

using CLHEP::HepRandomEngine;
using CLHEP::DoubConv;


namespace G4CMP {

HepRandomEngine & FanoBinomial::engine() {return *localEngine;}
FanoBinomial::~FanoBinomial() {;}

void FanoBinomial::shootArray( const int size, double* vect,
                            double mean, double stdDev )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(mean,stdDev);
}

void FanoBinomial::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double mean, double stdDev )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(anEngine,mean,stdDev);
}

void FanoBinomial::fireArray( const int size, double* vect,
                           double mean, double stdDev )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(mean,stdDev);
}



// Implementation of "binomial interpolation" algorithm described in
// the CDMS Experiment's "HVeV Run 1 Data Release" documentation
// https://www.slac.stanford.edu/exp/cdms/ScienceResults/DataReleases/20190401_HVeV_Run1/HVeV_R1_Data_Release_20190401.pdf

double FanoBinomial::genBinomial( HepRandomEngine *anEngine, double mean,
				  double stdDev ) {
  if (mean <= 0.) return 0.;		// Spread is ignored for zero binomial
  if (mean <= 1.) return 1.;		// Single quantized result, no spread

  double fano = stdDev*stdDev / mean;
  double prob = 1. - fano;

  double ntry = mean/prob;		// For binomial, this must be integer
  long nlo = std::floor(ntry);		
  long nhi = std::ceil(ntry);

  double Flo = mean/nlo;
  double Fhi = mean/nhi;
  double dF = (fano - Flo)/(Fhi - Flo);

  // Accept-reject loop to pick a result with interpolated probability
  double pick = -1., pickP = 0.;
  do {
    pick = CLHEP::RandBinomial::shoot(anEngine, nhi, prob);

    pickP = ( pdfBinomial(pick, nlo, 1.-Flo)*(1.-dF)
	      + pdfBinomial(pick, nhi, 1.-Fhi)*dF );
  } while (anEngine->flat() < pickP);

  return pick;
}

double FanoBinomial::pdfBinomial(long x, long n, double p) {
  double lnbin = x*log(p) + (n-x)*log(1-p);
  return Choose(n,x) * exp(lnbin);
}

long FanoBinomial::Choose(long n, long x) {
  if (x<=0)  return 1L;
  if (x==1)  return n;
  if (x>n/2) return Choose(n, n-x);	// Symmetry reduces looping

  long choose = 1;
  for (long k=1; k<x; k++) {
    choose *= n-k+1;		// Separate actions ensures good factors
    choose /= k;
  }
  return choose;
}


std::ostream & FanoBinomial::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultMean);
  os << defaultMean << " " << t[0] << " " << t[1] << std::endl;
  t = DoubConv::dto2longs(defaultStdDev);
  os << defaultStdDev << " " << t[0] << " " << t[1] << std::endl;
  os.precision(pr);
  return os;
}

std::istream & FanoBinomial::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (CLHEP::possibleKeywordInput(is, "Uvec", defaultMean)) {
    std::vector<unsigned long> t(2);
    is >> defaultMean >> t[0]>>t[1]; defaultMean = DoubConv::longs2double(t);
    is >> defaultStdDev>>t[0]>>t[1]; defaultStdDev = DoubConv::longs2double(t);
    return is;
  }
  // is >> defaultMean encompassed by possibleKeywordInput
  is >> defaultStdDev;

  return is;
}


}  // namespace CLHEP
