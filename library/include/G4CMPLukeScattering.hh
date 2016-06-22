/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLukeScattering.hh
/// \brief Definition of the G4CMPLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20160110  Remerge the electron and hole subclasses into one class

#ifndef G4CMPLukeScattering_h
#define G4CMPLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4ThreeVector.hh"
#include <iostream>

class G4CMPTrackInformation;
class G4VProcess;
class G4ParticleDefinition;
class G4Track;

class G4CMPLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPLukeScattering();
  virtual ~G4CMPLukeScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4double CalculateKSound(const G4CMPTrackInformation* trackInfo);
  G4double ComputeTimeSteps(const G4Track& aTrack);
  G4double TimeStepInField(G4double Emag, G4double coeff, G4double l0) const;
  G4double CalculateEmissionRate(const G4CMPTrackInformation* trackInfo, 
                                 const G4ThreeVector& k, 
                                 const G4ThreeVector& E, 
                                 G4double charge, 
                                 G4double t) const;
  // hide assignment operator as private
  G4CMPLukeScattering(G4CMPLukeScattering&);
  G4CMPLukeScattering& operator=(const G4CMPLukeScattering& right);

  std::vector<std::vector<std::vector<G4double> > > MeanFreeTimes;
  G4int ESIZE;
  G4int KSIZE;
  G4int THETASIZE;
  G4double EMAX;
  G4double KMAX;
  G4double THETAMAX;

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4CMPLukeScattering */
