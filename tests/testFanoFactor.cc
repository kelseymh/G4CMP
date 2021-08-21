// testFanoFactor: Exercise G4CMPFanoBinomial with non-integer N-true
//
// Usage: testFanoFactor Ntrue Fano
//
// Arguments: Ntrue can be a non-integer expected number of charge pairs,
//	      Fano is the Fano-factor desired.
//
// Code will generate 1000000 throws of G4CMPFanoBinomial, computing mean
// and standard deviation of output, to compare with input arguments.
//
// 20201213  Michael Kelsey
// 20210818  Add report of range of values

#include "globals.hh"
#include "G4CMPFanoBinomial.hh"
#include "Randomize.hh"
#include <float.h>
#include <stdlib.h>


// Throw 1M trials, get output mean and sigma
// NOTE: Using naive variance calculation, since use case is small Ntrue

void testFanoBinomial(double Ntrue, double Fano, G4int Ntrial) {
  G4double nthrow, mean=0., var=0.;
  G4double fbmin=DBL_MAX, fbmax=-DBL_MAX;
  for (G4int i=0; i<Ntrial; i++) {
    nthrow = G4CMP::FanoBinomial::shoot(Ntrue, Fano);
    mean += nthrow;
    var += nthrow*nthrow;

    if (nthrow < fbmin) fbmin = nthrow;
    if (nthrow > fbmax) fbmax = nthrow;
  }

  mean /= Ntrial;
  var = var/Ntrial - mean*mean;

  G4double fanoEst = var/mean;

  G4cout << "G4CMPFanoBinomial " << Ntrial << " throws:"
	 << " mean " << mean << " sigma " << sqrt(var) << " Fano " << fanoEst
	 << "\n range " << fbmin << " to " << fbmax
	 << G4endl;
}


// Copies of functions from G4CMPFanoBinomial, for validations

double Choose(long n, long x) {
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

double pdfBinomial(long x, long n, double p) {
  if (x<0 || x>n) return 0.;		// Non-physical
  if (p<=0. || p>1.) return 0.;
  if (p==1.) return (x==n)?1.:0.;	// Delta function at n

  double lnbin = log(p)*x + log(1.-p)*(n-x);	// p^x * (1-p)^(n-x)

  // Upper and lower limits for exp(), not defined anywhere else?
  const double DBL_EXPARG_MAX=709.7, DBL_EXPARG_MIN=-708.3;
  if (lnbin < DBL_EXPARG_MIN) return 0.;
  if (lnbin > DBL_EXPARG_MAX) return DBL_MAX;

  return Choose(n,x) * exp(lnbin);
}

double pdfInterp(double lo, double hi, double dP) {
  return (lo*(1.-dP) + hi*dP);
}

double binomInterp(long x, double mean, double prob) {
  double ntry = mean/prob;

  // Lower and upper bounds around non-integer "Ntry" value
  long nlo = std::floor(ntry);		
  long nhi = std::ceil(ntry);

  double Plo = mean/nlo;
  double Phi = mean/nhi;
  double dP = (prob - Plo)/(Phi - Plo);

  double pdfLo = pdfBinomial(x, nlo, Plo);
  double pdfHi = pdfBinomial(x, nhi, Phi);
  double pickP = pdfInterp(pdfLo, pdfHi, dP);

  return pickP;
}


// Validate FanoBinomial::pdfBinomial() calculation

void testPdfBinomial(long Ntrue, double Fano) {
  G4double prob = 1.-Fano;
  G4double Nmax = std::ceil(Ntrue/prob);

  G4double pdfMax = pdfBinomial(std::floor((Nmax+1)*prob), Nmax, prob);

  G4cout << "testPdfBinomial Nmax " << Nmax << " prob " << prob
	 << " pdfMax " << pdfMax << G4endl;

  for (long x=0; x<=Nmax; x++) {
    G4double pdf = pdfBinomial(x,Nmax,prob);
    G4cout << "testPdfBinomial " << x << " : " << pdf << " /pdfMax "
	   << pdf/pdfMax << G4endl;
  }
}


// Validate FanoBinomial interpolation expression for non-interger Ntry

void testInterp(double mean, double Fano) {
  long Nmax = std::ceil(mean/(1.-Fano));

  double pickP;
  double sumwt=0., sum=0., sumsq=0.;
  for (long x=0; x<=Nmax; x++) {
    pickP = binomInterp(x, mean, 1.-Fano);
    sumwt += pickP;		  // NOTE: Sum of weights should equal 1.
    sum   += x*pickP;
    sumsq += x*x*pickP;
  }

  sum /= sumwt;
  double var = sumsq/sumwt - sum*sum;
  double newFano = var/sum;

  G4cout << "testInterp mean " << sum << " stdDev " << sqrt(var)
	 << " Fano " << newFano << " sumwt " << sumwt << G4endl;
}


// Validate use of accept-reject loop

void testAcceptReject(double Ntrue, double Fano, G4int Ntrial) {
  std::map<long, long> samples;		// Key is x, value is count
  double prob = 1.-Fano;
  long Nmax = std::ceil(Ntrue/prob);

  double pdfMax = binomInterp(std::ceil(Ntrue), Ntrue, prob);
  G4cout << "testAcceptReject pdfMax " << pdfMax << G4endl;

  long pick;
  double pickP;
  G4int nthrow=0;
  for (G4int i=0; i<Ntrial; i++) {
    do {
      nthrow++;		// Count average throws required for a pick

      // Choose N value uniformly, or from binomial
      pick = std::floor((Nmax+1)*G4UniformRand());
      //*** pick = CLHEP::RandBinomial::shoot(Nmax, prob);

      pickP = binomInterp(pick, Ntrue, prob);

      if (pickP > pdfMax) {
	G4cerr << "WARNING! pick " << pick << " returns pickP " << pickP
	       << " > pdfMax " << pdfMax << G4endl;
      }
    } while (G4UniformRand() > pickP/pdfMax);

    samples[pick]++;		// Accumulate chosen values
  }

  double sum=0., sumsq=0.;
  for (long N=0; N<=Nmax; N++) {
    G4cout << "testAcceptReject " << N << " : " << samples[N] << " "
	   << double(samples[N])/double(Ntrial) << G4endl;
    sum += N*samples[N];
    sumsq += N*N*samples[N];
  }

  sum /= Ntrial;
  double var = sumsq/Ntrial - sum*sum;

  G4cout << "testAcceptReject " << double(nthrow)/Ntrial << " avg. throws:"
	 << " mean " << sum << " stdDev " << sqrt(var)
	 << " Fano " << var/sum << G4endl;
}


// MAIN PROGRAM

int main(int argc, char* argv[]) {
  // Get required arguments, or report usage and exit
  if (argc < 3) {
    G4cerr << "Usage: testFanoFactor Ntrue Fano [Nthrow]\n\n"
	   << "Arguments: Ntrue is expected number of charge pairs\n"
	   << "           Fano is Fano factor for material.\n"
	   << "           Nthrow is number of samples to generate.\n\n"
	   << "Code returns mean, sigma and computed Fano factor from"
	   << " Nthrow (10000) throws of\n G4CMPFanoBinomial." << G4endl;
    ::exit(1);
  }

  G4double Ntrue = strtod(argv[1],0);
  G4double Fano  = strtod(argv[2],0);
  G4int Nthrow = (argc>3) ? atoi(argv[3]) : 10000;

  G4cout << argv[0] << " " << Ntrue << " " << Fano << " " << Nthrow << G4endl;
  //*** testPdfBinomial(Ntrue, Fano);
  //*** testInterp(Ntrue, Fano);
  //*** testAcceptReject(Ntrue, Fano, Nthrow);
  testFanoBinomial(Ntrue, Fano, Nthrow);
}
