/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPSCUtils.cc
/// \brief Implementation of the G4CMPSCUtils class

#include "G4CMPSCUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"
#include "G4PhysicalConstants.hh"
#include "G4LatticePhysical.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor
// In principle, should run once per process. But in principle we can make static to reduce computation time?
// Here, since we're using a c-like array for efficiency, we have to initialize this in the constructor
G4CMPSCUtils::G4CMPSCUtils()
  :   fGapEnergyTempDependence { { 0.0, 1.0 },    
				 { 0.0015700754863653965, 0.9999846572428044 },
				 { 0.03768181167276968, 0.9996317738273021 },
				 { 0.07536362334553937, 0.999263547654604 },
				 { 0.12251191620808871, 0.9955454861608332 },
				 { 0.15705357690812766, 0.9952079455025266 },
				 { 0.1947814168524845, 0.9915819405519303 },
				 { 0.22618292657979266, 0.9912750854080152 },
				 { 0.257630464578688, 0.9877104514862018 },
				 { 0.2812276251457563, 0.9842225313503673 },
				 { 0.30796493668555547, 0.9807039257001413 },
				 { 0.32842194627989296, 0.9772466910786983 },
				 { 0.347308880387865, 0.973804799214451 },
				 { 0.3677658899822026, 0.970347564593008 },
				 { 0.4024456354970032, 0.9602366876010066 },
				 { 0.4324151545527076, 0.9501718388805924 },
				 { 0.4529181924186324, 0.9434568254812512 },
				 { 0.4813176359879714, 0.933407319518033 },
				 { 0.5034367776118489, 0.9234191845835977 },
				 { 0.5239858437493609, 0.913446392406358 },
				 { 0.5461049853732381, 0.9034582574719229 },
				 { 0.5714103062414339, 0.8901816582451978 },
				 { 0.5935294478653113, 0.8801935233107625 },
				 { 0.6236370517357775, 0.8603553382566538 },
				 { 0.6410919951721461, 0.8471554528159074 },
				 { 0.6633031933391978, 0.8306517603256757 },
				 { 0.679188061289201, 0.8174672176421252 },
				 { 0.6983051367551094, 0.7977364318883867 },
				 { 0.7158061084630651, 0.7812787676697421 },
				 { 0.7381093631732916, 0.7582595176237138 },
				 { 0.7557023914244219, 0.735286295849273 },
				 { 0.7716332876460119, 0.7188439743878241 },
				 { 0.7924124951414604, 0.6925822883210934 },
				 { 0.8100055233925907, 0.6696090665466523 },
				 { 0.8276906081868955, 0.6401202872164148 },
				 { 0.8469457684675656, 0.6106161651289815 },
				 { 0.8646768815334578, 0.5778696070208458 },
				 { 0.8808379191129843, 0.5451383916699057 },
				 { 0.8905345416607002,  0.5254996624593418 },
				 { 0.8986610887220507,  0.5058762760059735 },
				 { 0.913252050815212,  0.473160403412229 },
				 { 0.926272937422008, 0.4404598735756807 },
				 { 0.9394319088435656, 0.3979860074054375 },
				 { 0.9509747765071703, 0.3587852627702883 },
				 { 0.9642718327434898, 0.30653806026635033 },
				 { 0.9744747662786657, 0.25106376449890555 },
				 { 0.9813074074831742, 0.21190904813534361 },
				 { 0.9866620297444924, 0.16625411697318082 },
				 { 0.9920626802773973, 0.11734140703312002 },
				 { 0.9940930384796354, 0.08474827649694183 },
				 { 0.997739500439826, 0.04888202442566958 },
				 { 1.0, 0.0 } }

{
  G4cout << "REL HereA_SCUtils" << G4endl;
  
  //For now, (VERY TEMPORARY), define default SC parameters
  SetCurrentSCInfoToNull();

  //These are for setting ranges in our lookup tables. Not sure what the best bounds
  //are here yet, but tbd. Note: NEED TO SET THESE HIGHER THAN THE MAX EXPECTED PHONON ENERGY, OR RATES WILL BE WRONG.
  fMinPhononEnergyDivGap = 2.0;
  fMaxPhononEnergyDivGap = 60.0; //Was 20
  fMinQPEnergyDivGap = 1.0;
  fMaxQPEnergyDivGap = 60.0; //Was 20
  
  //Load "dimensionless" parameters that are universal to all BCS superconductors
  //  LoadGapEnergyTemperatureDependence();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Destructor
G4CMPSCUtils::~G4CMPSCUtils()
{  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Load this lattice's SC info into a local copy for convenience in calculations
void G4CMPSCUtils::LoadLatticeInfoIntoSCUtils(const G4LatticePhysical * theLat)
{
  fGap0Energy = theLat->GetSCDelta0();
  fTau0_qp = theLat->GetSCTau0qp();
  fTau0_ph = theLat->GetSCTau0ph();
  fTcrit = theLat->GetSCTcrit();
  fTeff = theLat->GetSCTeff();

  //REL. I think (?) we may want to move this one out to be in a dedicated rate process attached to the transport class.
  //It seems like it's somewhat specific to a single process. Ask Eric. I think what I want in here are parameters that are
  //useful for more than one class. This one seems a bit like the polycrystal elastic scattering mfp, which doesn't need to
  //be updated here because it's handled in the lattice itself and isn't part of any intensive calculations that need doing.
  fDn = theLat->GetSCDn();    
  
  fGapEnergy = ComputeCurrentGapEnergyAtNonzeroT();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Make parameters null so that we don't accidentally run any processes in silicon
void G4CMPSCUtils::SetCurrentSCInfoToNull()
{
  fGap0Energy = 0;
  fTau0_qp = DBL_MAX;
  fTau0_ph = DBL_MAX;
  fTcrit = 0;
  fTeff = 0;
  fDn = 0;
  fGapEnergy = 0;  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Take the gap energy at T=0, combine with knowledge of Teff/Tcrit, and compute a temperature-dependent gap.
G4double G4CMPSCUtils::ComputeCurrentGapEnergyAtNonzeroT()
{
  //Here, we use a simplified assumption, which is that the gap is that corresponding to T=Teff.
  double TeffDivTcrit = fTeff/fTcrit;
  if( TeffDivTcrit >=1 ){
    G4ExceptionDescription msg;
    msg << "Attempting to compute gap energy at Teff=" << fTeff << " while Tcrit=" << fTcrit << ". Whatever you think is superconducting is now not superconducting... Something might be wrong!";
    G4Exception("G4CMPSCUtils::ComputeCurrentGapEnergyAtNonzeroT", "G4CMPSCUtils001",JustWarning, msg);
    return 0;
  }
  else{
    //Now actually do the lookup. Assuming we're only establishing a temperature once per instantiation of a new SC, we only have to do this once.
    //To make it straightforward, we can do a loop (but this is a place where we can refine for speed, REL)
    double gapFactor = 0;
    for( int iT = 0; iT < fGapEnergyTempDependenceBins-1; ++iT ){
      if( (fTeff/fTcrit) >= fGapEnergyTempDependence[iT][0] && (fTeff/fTcrit) < fGapEnergyTempDependence[iT+1][0] ){
	gapFactor = fGapEnergyTempDependence[iT][1];
      }
    }
    return (fGap0Energy*gapFactor);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//In some cases, we actually want to be able to compute the gap at a nonzero T for volumes that the particle undergoing
//a process (which derives from this class) hasn't entered yet, but needs to check. (Ex: QPs checking gap of nearby volume to
//determine if transmission occurs). Since that Teff/Tcrit info hasn't been loaded into the internal data members here,
//we want a version of the above function that can be callable externally given an "external" gap, Teff, and Tcrit
G4double G4CMPSCUtils::ComputeTestGapEnergyAtNonzeroT(double Teff, double Tcrit, double gap0Energy) const
{
  //Here, we use a simplified assumption, which is that the gap is that corresponding to T=Teff.
  double TeffDivTcrit = Teff/Tcrit;
  if( TeffDivTcrit >=1 ){
    G4ExceptionDescription msg;
    msg << "Attempting to compute test gap energy at Teff=" << Teff << " while Tcrit=" << Tcrit << ". Whatever you think is superconducting is now not superconducting... Something might be wrong!";
    G4Exception("G4CMPSCUtils::ComputeTestGapEnergyAtNonzeroT", "G4CMPSCUtils002",JustWarning, msg);
    return 0;
  }
  else{
    //Now actually do the lookup. Assuming we're only establishing a temperature once per instantiation of a new SC, we only have to do this once.
    //To make it straightforward, we can do a loop (but this is a place where we can refine for speed, REL). We also can use the internal
    //fGapEnergyTempDependence because it's normalized, and is instantiated the same regardless of what the fTeff and fTcrit are
    double gapFactor = 0;
    for( int iT = 0; iT < fGapEnergyTempDependenceBins-1; ++iT ){
      if( (Teff/Tcrit) >= fGapEnergyTempDependence[iT][0] && (Teff/Tcrit) < fGapEnergyTempDependence[iT+1][0] ){
	gapFactor = fGapEnergyTempDependence[iT][1];
      }
    }
    return (gap0Energy*gapFactor);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4CMPSCUtils::FermiFactor(G4double energy, G4double temperature)
{
  return 1/(exp((energy)/(k_Boltzmann * temperature))+1);  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4CMPSCUtils::BoseFactor(G4double energy, G4double temperature)
{
  return 1/(exp((energy)/(k_Boltzmann * temperature))-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Get a Tau from a lookup table vs. energy. Note that the table has two columns: energy and tau. For now we ignore the former but we
// may end up using this to help check ourselves. Passing by const reference
G4double G4CMPSCUtils::GetTauAsAFunctionOfEnergy( const std::vector<std::vector<G4double> > & tauVsPhononEnergy, G4String particleInQuestion, G4double energy, G4bool & thisEnergyBelowUsableRange ) const
{
  //Establish what the bounds and binning are
  G4double minE = 0;
  G4double maxE = 0;
  G4int nE = 0;
  if( particleInQuestion == "Phonon" ){
    minE = fMinPhononEnergyDivGap * fGapEnergy;
    maxE = fMaxPhononEnergyDivGap * fGapEnergy;
    nE = fPhononEnergyBins;
  }
  else if( particleInQuestion == "BogoliubovQP" ){
    minE = fMinQPEnergyDivGap * fGapEnergy;
    maxE = fMaxQPEnergyDivGap * fGapEnergy;
    nE = fQPEnergyBins;
  }
  else{
    G4ExceptionDescription msg;
    msg << "During lookup table step, particle in question is neither a phonon nor a BogoliubovQP. This is probably a bug in the G4CMPSCUtils code somewhere.";
    G4Exception("G4CMPSCUtils::GetTauAsAFunctionOfEnergy", "G4CMPSCUtils003",FatalException, msg);
    return 0;
  }    

  //Now we can guess at a location within the lookup table. First, handle some edge cases.
  if( energy <= minE ){
    
    //Split this into two edge cases: if we're looking at phonons, then it's okay for the energy to be below the min -- the lookup table should just
    //return zero.
    if( particleInQuestion == "Phonon" ){
      return DBL_MAX;
    }
    else if( particleInQuestion == "BogoliubovQP" ){

      //Should split this again here into "less than" and "equal to" -- the latter is in principle not physically inaccurate, but it seems
      //like floating point errors may cause it to occur? In any case, I'll allow this one as a test and return the tau at the lowest energy
      //bin. (REL take a look at this later -- it seems like this happens but the chance of a floating-point error to below the gap is
      //also seemingly negligible. Weird.)
      if( energy == minE ){
	return tauVsPhononEnergy[0][1];	
      }
      else{
	if( G4CMPConfigManager::GetVerboseLevel() > 5 ){
	  G4cout << "energy: " << energy << ", minE: " << minE << G4endl;
	  G4ExceptionDescription msg;
	  msg << "During lookup table step, we're somehow looking at a QP energy below the gap? This is probably just a turnaround step.";
	  G4Exception("G4CMPSCUtils::GetTauAsAFunctionOfEnergy", "G4CMPSCUtils004",JustWarning, msg);
	}
	thisEnergyBelowUsableRange = true;
	return DBL_MAX;
      }
    }
    else{
      return 0; //Shouldn't get here.
    }
  }

  //If we're above the maximum energy considered, for now let's just set it to the max bin value (which should basically be zero if the
  //bounds are set correctly). We should probably refine this a bit (REL) but that refinement may depend on which kind of particle AND
  //which kind of process we're looking at
  else if( energy >= maxE ){
    return tauVsPhononEnergy[nE-1][1];
  }

  //Now we do the guess. Note that this currently doesn't interpolate -- it just selects a single bin. REL Using the additional info
  //provided by the energy column, we should interpolate for best results.
  //As a good reminder, the min and max are taken to be the "bookends" to the energy range, but the actual bin/array values are
  //associated with energies that are offset by 0.5*deltaE (where deltaE is whatever width goes into that energy range corresponding
  //to the required binning). This is done to prevent any energy bins being exactly at the gap. (This maaaay not be 100% necessary,
  //but I'll need to rethink about it.)
  else{
    double deltaE = tauVsPhononEnergy[1][0]-tauVsPhononEnergy[0][0];
    int binGuess = floor((energy-minE) / deltaE);

    //If our bin guess is the final bin, then just use that bin value. We'll accept a single edge case here.
    if( binGuess == tauVsPhononEnergy.size()-1 ){
      return tauVsPhononEnergy[binGuess][1];
    }
    //Otherwise, try to interpolate.
    else{
    
      //Can compute the slope using the energy column and nearby bin values. Compute deltaE from these 
      double lowerE = tauVsPhononEnergy[binGuess][0];
      double upperE = tauVsPhononEnergy[binGuess+1][0];
      double lowerTau = tauVsPhononEnergy[binGuess][1];
      double upperTau = tauVsPhononEnergy[binGuess+1][1];

      
      //Cross-check. Noting that here, again, the E values of the array bins are set in the *middle* of the bins. So since our bin guess
      //uses a floor mechanism considering the whole range and deltaE, we have to shift down in energy by 0.5 deltaE for this
      if( energy >= upperE-0.5*deltaE || energy < lowerE-0.5*deltaE ){
	G4ExceptionDescription msg;
	msg << "During interpolation in lookup table step, we're somehow off-by-one in our energy bin... A calculation is being done wrong. Investigate. Energy: " << energy << ", UpperE: " << upperE << ", lowerE: " << lowerE;
	G4Exception("G4CMPSCUtils::GetTauAsAFunctionOfEnergy", "G4CMPSCUtils005",FatalException, msg);
      }
      
      double slope = (upperTau-lowerTau)/(upperE-lowerE);
      double interpolatedVal = lowerTau + slope*(energy-lowerE);
      return interpolatedVal;
    }
  }
}








//------------------------ BELOW THIS LINE NEED TO BE PLACED INTO CLASSES THAT HAVEN'T BEEN CONSTRUCTED YET SO ARE HERE FOR STORAGE
//--------------------------------------------------------------------------------------------------------------------------------------


/*
// DOES NOT REMAIN -- GOES INTO THE DEDICATED QP SCATTERING CLASS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for QP radiation vs. QP energy
G4double * G4CMPSCUtils::ComputeNormalizedTauQPRadiationVsEnergy()
{
  G4double output[fQPEnergyBins] = {0};
  G4double deltaQPEnergyDivGap = (fMaxQPEnergyDivGap - fMinQPEnergyDivGap) / ((double)fQPEnergyBins);

  //Loop over all QP energy bins, and create a normalized tau for each of them
  for( int iB = 0; iB < fQPEnergyBins; ++iB ){

    //Define the qp energy for this bin    
    G4double qpEnergy = fGapEnergy * (deltaQPEnergyDivGap * iB + fMinQPEnergyDivGap);
    
    //--------------------------------------------------
    //Now do the calculation for this phonon energy       
    //Initially, going to try this without specifying units -- shouuuuuld be fine but if there are problems this is the first place to debug.
    //Note that in what follows, Omega indicates the phonon energy integrated over.
    double prefactor = 2 * pi / hbar_Planck / (1-FermiFactor(qpEnergy,fTeff)); //Noting that there's a factor of Z1_0 missing here but there's a similar factor in tau0 that we're omitting as they cancel.
    int nW = 10000;
    double minOmega = 0;
    double maxOmega = qpEnergy - fGapEnergy;
    double deltaOmega = (maxOmega-minOmega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++ iW ){
      double Omega = minOmega + (iW+0.5)*deltaOmega;
      double alpha2F = Omega * Omega; //Noting that there's a factor of b missing here, but there is a 1/b in tau0 that we're also going to omit as they cancel.
      double energyTerm1 = (qpEnergy - Omega) / pow( pow(qpEnergy-Omega,2) - fGapEnergy*fGapEnergy, 0.5);
      double energyTerm2 = (1 - fGapEnergy*fGapEnergy / qpEnergy / (qpEnergy - Omega));
      double integrand = alpha2F * energyTerm1 * energyTerm2 * (BoseFactor(Omega,fTeff) + 1) * (1 - FermiFactor(qpEnergy-Omega,fTeff));
      integral += integrand * deltaOmega;
    }

    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 2 / pi / pow(k_Boltzmann*fTcrit,3); //Here, omitting Z1 and b a la the above comments
    double normalizedTau = 1 / inverseTau / tau0;

    output[iB] = { qpEnergy, normalizedTau };
  }
  return output;
}

// DOES NOT REMAIN -- GOES INTO THE DEDICATED QP SCATTERING CLASS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for QP-phonon scattering
G4double * G4CMPSCUtils::ComputeNormalizedTauQPScattersOnPhononVsEnergy()
{
  G4double output[fQPEnergyBins] = {0};
  G4double deltaQPEnergyDivGap = (fMaxQPEnergyDivGap - fMinQPEnergyDivGap) / ((double)fQPEnergyBins);
    
  //Loop over all QP energy bins, and create a normalized tau for each of them
  for( int iB = 0; iB < fQPEnergyBins; ++iB ){

    //Define the qp energy for this bin    
    G4double qpEnergy = fGapEnergy * (deltaQPEnergyDivGap * iB + fMinQPEnergyDivGap);
    
    //--------------------------------------------------
    //Now do the calculation for this phonon energy       
    //Initially, going to try this without specifying units -- shouuuuuld be fine but if there are problems this is the first place to debug.
    //Note that in what follows, Omega indicates the phonon energy integrated over.  
    double prefactor = 2 * pi / hbar_Planck / (1-FermiFactor(qpEnergy,fTeff)); //Omitting Z1_0 here b/c it's canceled in denom
    
    //Now do an integral over Omega, phonon energy
    int nW = 10000;
    double minOmega = 0; //Start with actually zero phonon energy
    double maxOmega = 0.1 * CLHEP::eV; //Go up to 100 meV in phonon energy -- really shouldn't need to go higher than this.
    double deltaOmega = (maxOmega - minOmega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++iW ){
      double Omega = minOmega + (iW+0.5)*deltaOmega;
      double alpha2F = Omega * Omega; //Noting that there's a factor of b missing here, but there is a 1/b in tau0 that we're also going to omit as they cancel.
      double energyTerm1 = (qpEnergy+Omega)/pow(pow(qpEnergy+Omega,2)-fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 = (1 - fGapEnergy*fGapEnergy/qpEnergy/(qpEnergy+Omega));
      double integrand = alpha2F * energyTerm1 * energyTerm2 * (BoseFactor(Omega,fTeff)) * (1 - FermiFactor(qpEnergy+Omega,fTeff));
      integral += integrand * deltaOmega;
      
      //std::cout << "PhononScatter integrand: " << integrand << std::endl;
    }    
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 2 / pi / pow(k_Boltzmann*fTcrit,3);  //Here, omitting Z1 and b a la the above comments
    double normalizedTau = 1/inverseTau / tau0;  
    output[iB] = { qpEnergy, normalizedTau };
  }
  return output;  
}


*/













/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Compute phonon energy distribution from quasiparticle in superconductor.
// This is the same as the original KaplanQP
G4double G4CMPSCUtils::PhononEnergyRand(G4double Energy) const
{
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  //
  // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-fGapEnergy*fGapEnergy/Energy))
  //           /
  //           sqrt((E'*E' - fGapEnergy*fGapEnergy);
  
  // Add buffer so first bin doesn't give zero denominator in pdfSum
  
  const G4double BUFF = 10000000000.; //REL used to be 1000 by default
  G4double xmin = fGapEnergy + fGapEnergy / BUFF;
  G4double xmax = Energy;
  G4double ymax = PhononEnergyPDF(Energy, xmin);
  
  G4double xtest = 0., ytest = ymax;
  do {
    ytest = G4UniformRand() * ymax;
    xtest = G4UniformRand() * (xmax - xmin) + xmin;
  } while (ytest > PhononEnergyPDF(Energy, xtest));  
  return Energy - xtest;
}

// DOES NOT REMAIN -- GOES INTO THE DEDICATED QP SCATTERING CLASS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Phonon energy "pdf"
// This is the same as the original KaplanQP, but with a modification to fix a typo and
// bring this into agreement with the Kaplan paper
G4double G4CMPSCUtils::PhononEnergyPDF(G4double E, G4double x) const
{
    const G4double gapsq = fGapEnergy * fGapEnergy;
    return (x * (E - x) * (E - x) * (1 - gapsq / x / E) / sqrt(x * x - gapsq));
}
*/

























//------------------------ BELOW THIS LINE SHOULD GET TRASHED ONCE WE CONFIRM THINGS ARE WORKING
//------------------------------------------------------------------------------------------------------------






/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Construct the lookup table for normalized tau for pairbreaking vs phonon energy
G4double * G4CMPSCUtils::ComputeNormalizedTauPairBreakingVsEnergy()
{
  G4double output[fPhononEnergyBins] = {0};
  G4double deltaPhononEnergyDivGap = (fMaxPhononEnergyDivGap - fMinPhononEnergyDivGap) / ((double)fPhononEnergyBins);  
  
  //Loop over all phonon energy bins, and for each of them, compute the integral
  for( int iB = 0; iB < fPhononEnergyBins; ++iB ){

    //Define the phonon energy for this bin    
    G4double phononEnergy = fGapEnergy * (deltaPhononEnergyDivGap * iB + fMinPhononEnergyDivGap);

    //--------------------------------------------------
    //Now do the calculation for this phonon energy       
    //Initially, going to try this without specifying units -- shouuuuuld be fine but if there are problems this is the first place to debug.
    //Note that in what follows, omega indicates the QP energy integrated over.
    double prefactor = 4 * pi / hbar_Planck * fAlpha2;
    int nW = 10000;
    double minomega = fGapEnergy;
    double maxomega = phononEnergy - fGapEnergy;
    double deltaomega = (maxomega - minomega) / ((double)nW);
    double integral = 0;
    for( int iW = 0; iW < nW; ++iW ){
      double omega = minomega + (iW+0.5)*deltaomega; //So that omega lands in the middle of a bin
      double energyTerm1 = 1.0 / pow(omega*omega - fGapEnergy*fGapEnergy,0.5);
      double energyTerm2 = (omega*(phononEnergy-omega) + fGapEnergy*fGapEnergy) / pow( pow(phononEnergy-omega,2) - fGapEnergy*fGapEnergy, 0.5 );
      double integrand = energyTerm1 * energyTerm2 * (1-FermiFactor(phononEnergy-omega,fTeff) - FermiFactor(omega,fTeff));
      integral += integrand * deltaomega;
    }
    double inverseTau = prefactor * integral;
    double tau0 = hbar_Planck / 4 / pi / pi / fAlpha2 / fGap0Energy;
    double normalizedTau = 1.0 / inverseTau / tau0;

    output[iB] = { phononEnergy, normalizedTau };
  }
  return output;
}
*/


















/*
//DOES NOT REMAIN
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Use the current lattice to update the SC utils class in various ways:
void G4CMPSCUtils::UpdateSCUtilsFunctions(const G4LatticePhysical* theLat)
{
  //First check: if the lattice doesn't represent a superconductor, then we need to
  //ensure that the parameters are set to "null" values and return -- don't want to do
  //any additional updates.
  if( !theLat->GetSCDelta0() <= 0.0 ){
    SetCurrentSCInfoToNull();
    return;
  }
  
  //1. If the lattice doesn't exist in the lattice container associated with this process yet,
  //   add it and do the full calculation of the curves we care about, storing them in a map
  if( fMap_physicalLattice_alreadyBuilt.count(theLat) == 0 ){
    fMap_physicalLattice_alreadyBuilt.emplace(theLat,true);
    LoadLatticeInfoIntoSCUtils(theLat);
    GenerateAllLifetimesVsEnergy(theLat);
    UpdateAllLifetimesVsEnergy(theLat);  
  }  
  
  //2. If it does, just make the "active" functions the ones in the map
  else{
    LoadLatticeInfoIntoSCUtils(theLat);
    UpdateAllLifetimesVsEnergy(theLat);  
  }
}
*/




/*
// DOES NOT REMAIN -- SPLIT INTO THE VARIOUS RATE CLASSES
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Generate the lifetime arrays/lookup tables for all of the various processes
// for this polycrystal. In principle there should never be a time when only a subset
// of these are generated.
void G4CMPSCUtils::GenerateAllLifetimesVsEnergy(G4LatticePhysical* theLat)
{
  //The computations here take the "current" SC parameters, which are set for this SC prior to running these functions.
  fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy.emplace(theLat,ComputeNormalizedTauPairBreakingVsEnergy());
  fMap_physicalLattice_NormalizedTauQPRadiationVsEnergy.emplace(theLat,,ComputeNormalizedTauQPRadiationVsEnergy());
  fMap_physicalLattice_NormalizedTauQPScattersOnPhononVsEnergy.emplace(theLat,ComputeNormalizedTauQPScattersOnPhononVsEnergy());
  fMap_physicalLattice_NormalizedTauQPRecombinationVsEnergy.emplace(theLat,ComputeNormalizedTauQPRecombinationVsEnergy()); 
}

// DOES NOT REMAIN -- SPLIT INTO THE VARIOUS RATE CLASSES
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Make sure to update the lifetime arrays/lookup tables for this polycrystal.
void G4CMPSCUtils::UpdateAllLifetimesVsEnergy(G4LatticePhysical* theLat)
{
  fNormalizedTauPairBreakingVsEnergy = fMap_physicalLattice_NormalizedTauPairBreakingVsEnergy[theLat];
  fNormalizedTauQPRadiationVsEnergy = fMap_physicalLattice_NormalizedTauQPRadiationVsEnergy[theLat];  
  fNormalizedTauQPScattersOnPhononVsEnergy = fMap_physicalLattice_NormalizedTauQPScattersOnPhononVsEnergy[theLat];
  fNormalizedTauQPRecombinationVsEnergy = fMap_physicalLattice_NormalizedTauQPRecombinationVsEnergy[theLat];
}



// DOES NOT REMAIN -- GOES INTO THE DEDICATED PAIRBREAKING CLASS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//For a given phonon energy, sample QP energies from the relevant distribution. 
std::pair<G4double,G4double> G4CMPSCUtils::FetchQPEnergies(G4double phononEnergy)
{
  //From the G4CMPSUtils class, sample a qp energy:
  G4double qp1Energy = QPEnergyRand(energy);
  G4double qp2Energy = phonEnergy - qp1Energy;
  std::pair<G4double,G4double> thePair(qp1Energy,qp2Energy);
  return thePair;
  
}


*/

