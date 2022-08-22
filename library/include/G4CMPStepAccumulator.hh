/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPStepAccumulator.hh
/// \brief Definition of the G4CMPStepAccumulator container.  This class  
///	stores energy deposit information from track steps, to create
///     phonon and charge carrier secondaries.
//
// 20210303  Michael Kelsey
// 20210608  Ensure that all values are initialized in constructor; add
//	       eventID to avoid rollover between events.
// 20220216  Add "stepID" to provide full identification.  Add printout.
// 20220821  G4CMP-308 -- Define step-info container to avoid needing G4Step

#ifndef G4CMPStepAccumulator_hh
#define G4CMPStepAccumulator_hh 1

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>

class G4ParticleDefinition;
class G4Step;


class G4CMPStepInfo {
public:
  G4CMPStepInfo() : trackID(-1), stepID(-1), pd(0), Edep(0.), Eniel(0.) {;}
  G4CMPStepInfo(const G4Step* step);
  G4CMPStepInfo(const G4Step& step);
  G4CMPStepInfo& operator=(const G4CMPStepInfo& step);

  virtual void Clear() {
    trackID = stepID = -1;
    pd = 0; Edep = Eniel = 0.;
    start.set(0,0,0); end.set(0,0,0);
  }

public:
  G4int trackID;		  // Track ID for sanity checking in Add()
  G4int stepID;		  	  // Step ID of last step accumulated
  const G4ParticleDefinition* pd; // Particle type (should match all steps)
  G4double Edep;		  // Sum of GetTotalEnergyDeposit()
  G4double Eniel;		  // Sum of GetNonIonizingEnergyDeposit()
  G4ThreeVector start;	  	  // PreStep position of first hit
  G4ThreeVector end;		  // PostStep position of last hit

};

class G4CMPStepAccumulator : public G4CMPStepInfo {
public:
  G4CMPStepAccumulator() : G4CMPStepInfo(), nsteps(0), eventID(-1) {;}
  ~G4CMPStepAccumulator() {;}

  // Extract relevant information from step
  void Add(const G4CMPStepInfo& step);
  void Add(const G4Step& step) { Add(G4CMPStepInfo(step)); }
  void Add(const G4Step* step) { if (step) Add(*step); }

  // Reset accumulator for new track
  virtual void Clear() { G4CMPStepInfo::Clear(); nsteps=0; eventID=-1; }

  // Dump content for diagnostics
  void Print(std::ostream& os) const;

public:
  G4int nsteps;			// Accumulated steps so far
  G4int eventID;		// Current event ID for sanity checking in Add()
};

inline 
std::ostream& operator<<(std::ostream& os, const G4CMPStepAccumulator& accum) {
  accum.Print(os);
  return os;
}

#endif	/* G4CMPStepAccumulator_hh */
