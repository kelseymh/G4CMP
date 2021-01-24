// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                        --- FanoBinomial ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// 20201019  Michael Kelsey (TAMU)
// 20201213  Replace "stdDev" argument with Fano factor
// 20201218  Improve computing steps to avoid loss of precision, apply
//	       cutoff for binomial -> Gaussian limit.
// 20210107  CLHEP's DoubConv.* header file has different names in CLHEP
//		distribution vs. Geant4 internal subset.
// 20210123  Strip all use of DoubConv (broken for us in CLHEP 2.4.4.1)
// =======================================================================

#include "G4CMPFanoBinomial.hh"
#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include <algorithm>	// for min() and max()
#include <cfloat>	// for DBL_MAX
#include <cmath>	// for exp()
#include <iostream>

using CLHEP::HepRandomEngine;


namespace G4CMP {

HepRandomEngine & FanoBinomial::engine() {return *localEngine;}
FanoBinomial::~FanoBinomial() {;}

void FanoBinomial::shootArray( const int size, double* vect,
                            double mean, double fano )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(mean,fano);
}

void FanoBinomial::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double mean, double fano )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(anEngine,mean,fano);
}

void FanoBinomial::fireArray( const int size, double* vect,
                           double mean, double fano )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(mean,fano);
}



// Implementation of "binomial interpolation" algorithm described in
// the CDMS Experiment's "HVeV Run 1 Data Release" documentation
// https://www.slac.stanford.edu/exp/cdms/ScienceResults/DataReleases/20190401_HVeV_Run1/HVeV_R1_Data_Release_20190401.pdf

double FanoBinomial::genBinomial( HepRandomEngine *anEngine, double mean,
				  double fano ) {
  if (mean <= 0.) return 0.;		// Spread is ignored for zero binomial
  if (mean <= 1.) return 1.;		// Single quantized result, no spread

  // Fano factor of 1. means Poisson distribution
  if (fano == 1.) return CLHEP::RandPoissonQ::shoot(anEngine, mean);

  double prob = 1. - fano;
  double ntry = mean/prob;

  // Use Gaussian approximation where appropriate
  if (mean > 9*(fano/prob) && mean > 9*(prob/fano))
    return CLHEP::RandGaussQ::shoot(anEngine, mean, sqrt(mean*fano));

  // If integer, then nlo==nhi and no interpolation is needed
  if (ntry == int(ntry)) {
    return CLHEP::RandBinomial::shoot(anEngine, long(ntry), prob);
  }

  // Implement interpolated binomial distribution
  long nlo = std::floor(ntry);		
  long nhi = std::ceil(ntry);

  double Plo = mean/nlo;
  double Phi = mean/nhi;
  double dP = (prob - Plo)/(Phi - Plo);

  // Accept-reject loop to pick a result with interpolated probability
  double pdfMax = pdfBinomial(long(nhi*prob), nlo, Plo);

  long pick;
  double pickP;
  do {
    pick = long(anEngine->flat()*(nhi+1));
    //*** pick = CLHEP::RandBinomial::shoot(anEngine, nhi, prob);

    pickP = ( pdfBinomial(pick, nlo, Plo)*(1.-dP)
	      + pdfBinomial(pick, nhi, Phi)*dP );
    pickP /= pdfMax;
  } while (anEngine->flat() > pickP);

  return pick;
}

double FanoBinomial::pdfBinomial(long x, long n, double p) {
  if (x > n || x < 0)  return 0.;		// Non-physical arguments
  if (p <= 0. || p>1.) return 0.;
  if (p == 1.) return (x==n ? 1. : 0.);		// Delta function at n

  double lnbin = log(p)*x + log(1.-p)*(n-x);	// p^x * (1-p)^(n-x)

  // Upper and lower limits for exp(), not defined anywhere else?
  const double DBL_EXPARG_MAX=709.7, DBL_EXPARG_MIN=-708.3;
  if (lnbin < DBL_EXPARG_MIN) return 0.;
  if (lnbin > DBL_EXPARG_MAX) return DBL_MAX;

  return Choose(n,x) * exp(lnbin);
}

double FanoBinomial::Choose(long n, long x) {
  if (x>n/2) x = n-x;			// Symmetry reduces looping
  if (x<0)  return 0.;			// Non-physical
  if (x==0) return 1.;			// Simple cases avoid computation
  if (x==1) return n;

  double choose = 1.;
  for (long k=1L; k<=x; k++) {
    choose *= double(n-k+1L)/k;
  }
  return choose;
}


std::ostream & FanoBinomial::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  os << defaultMean << std::endl;
  os << defaultFano << std::endl;
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
    is >> defaultMean >> defaultFano;
    return is;
  }
  // is >> defaultMean encompassed by possibleKeywordInput
  is >> defaultFano;

  return is;
}


}  // namespace CLHEP
