/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPInterValleyRate.hh
/// \brief Compute electron intervalley scattering rate using matrix elements
//
// $Id$
//
// 20170919  Add interface for threshold identification

#ifndef G4CMPInterValleyRate_hh
#define G4CMPInterValleyRate_hh 1

#include "G4CMPVScatteringRate.hh"


class G4CMPInterValleyRate : public G4CMPVScatteringRate {
public:
  G4CMPInterValleyRate()
    : G4CMPVScatteringRate("InterValley"),
      hbar_sq(CLHEP::hbar_Planck*CLHEP::hbar_Planck), hbar_4th(hbar_sq*hbar_sq),
      m_electron(CLHEP::electron_mass_c2/CLHEP::c_squared),
      eTrk(0.), density(0.), uSound(0.), alpha(0.),
      m_DOS(0.), m_DOS3half(0.) {;}

  virtual ~G4CMPInterValleyRate() {;}

  virtual G4double Rate(const G4Track& aTrack) const;

  virtual G4double Threshold(G4double Eabove=0.) const;

  // Initialize numerical parameters below
  virtual void LoadDataForTrack(const G4Track* track);
    
  //G4CMPInterValleyRate& operator=(const G4CMPInterValleyRate& rhs2);  
  const std::vector<G4double>& GetIVProb() const {return IVprob;} 
  
    

protected:
  G4double opticalRate() const;		// Optical intervalley D0, D1 rate
  G4double scatterRate() const;		// Neutral impurity scattering

  G4double energyFunc(G4double E) const {	// Energy dependence of rates
    return sqrt(E*(1+alpha*E))*(1+2*alpha*E);
  }

private:
  // Useful numerical parameters for computing individual rates
  const G4double hbar_sq;
  const G4double hbar_4th;
  const G4double m_electron;

  // Kinematic parameters set by Rate()
  mutable G4double eTrk;	// Track kinetic energy
  mutable G4int ivalley;	// Valley's index
  mutable G4ThreeVector ptrk;	// Local momentum
  mutable G4ThreeVector ktrk;	// Local quasi-momentum
  mutable G4ThreeVector kHV;	// quasi-momentum in HV frame
  mutable G4double kmag;	// magnitude of kHV
  mutable std::vector<G4double> IVprob; 	// Store IV rates

  G4double density;		// Crystal density (from G4Material)
  G4double uSound;		// Average sound speed for acoustic rate
  G4double alpha;		// Non-parabolicity parameter
  G4double m_DOS;		// Electron "density of states" average mass
  G4double m_DOS3half;		// m_DOS ^ (3/2)
  
};

#endif	/* G4CMPInterValleyRate_hh */
