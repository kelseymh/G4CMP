/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPFanoBinomial.hh
/// \brief Random number generator for a binomial distribution represented by
///	   a mean and a standard deviation (Fano factor).  Interpolates
///	   between adjacent integer 'N' binomial distributions with associated
///	   'p' probabilities derived from the input mean and sigma.
//
// 20201018  Michael Kelsey (TAMU) 

#ifndef G4CMPFanoBinomial_h
#define G4CMPFanoBinomial_h 1

#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"
#include <string>


namespace G4CMP {

class FanoBinomial : public CLHEP::HepRandom {
public:
  inline FanoBinomial ( CLHEP::HepRandomEngine& anEngine, double nTrue=0.,
			double fano=1. );
  inline FanoBinomial ( CLHEP::HepRandomEngine* anEngine, double nTrue=0.,
			double fano=1. );
  // These constructors should be used to instantiate a Fano-factor binomial
  // distribution object defining a local engine for it.  The Fano factor
  // F = stdDev^2 / mean.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the FanoBinomial destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the FanoBinomial destructor.

  virtual ~FanoBinomial();
  // Destructor

  // Static methods to shoot random values using the static generator

  static  inline double shoot( double mean, double fano );

  static  void shootArray ( const int size, double* vect,
                            double mean=0.0, double fano=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  double shoot( CLHEP::HepRandomEngine* anEngine );

  static  inline double shoot( CLHEP::HepRandomEngine* anEngine, 
			       double mean, double fano );

  static  void shootArray ( CLHEP::HepRandomEngine* anEngine, const int size,
                            double* vect, double mean=0.0,
                            double fano=1.0 );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline double fire();

  double fire( double mean, double fano );
  
  void fireArray ( const int size, double* vect);
  void fireArray ( const int size, double* vect,
                   double mean, double fano );

  inline double operator()();
  inline double operator()( double mean, double fano );

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  std::string name() const { return "FanoBinomial"; }
  CLHEP::HepRandomEngine & engine();

  static std::string distributionName() { return "FanoBinomial"; }  
  // Provides the name of this distribution class

private:
  static double genBinomial( CLHEP::HepRandomEngine *anEngine,
			     double mean, double fano );

  static double pdfBinomial(long x, long n, double p);

  static double Choose(long n, long x);

  std::shared_ptr<CLHEP::HepRandomEngine> localEngine;
  double defaultMean;
  double defaultFano;
};

}	// namespace G4CMP

#include "G4CMPFanoBinomial.icc"

#endif	/* G4CMPFanoBinomial_h */
