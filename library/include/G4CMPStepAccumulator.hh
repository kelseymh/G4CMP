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

#ifndef G4CMPStepAccumulator_hh
#define G4CMPStepAccumulator_hh 1

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>

class G4ParticleDefinition;
class G4Step;


class G4CMPStepAccumulator {
public:
  G4CMPStepAccumulator() { Clear(); }
  ~G4CMPStepAccumulator() {;}

  // Extract relevant information from step
  void Add(const G4Step& step);
  void Add(const G4Step* step) { if (step) Add(*step); }

  // Reset accumulator for new track
  void Clear();

  // Dump content for diagnostics
  void Print(std::ostream& os) const;

public:		// Data values are public for more convenient access
  G4int eventID;		  // Event ID for sanity checking in Add()
  G4int trackID;		  // Track ID for sanity checking in Add()
  G4int stepID;			  // Step ID of last step accumulated
  G4int nsteps;			  // Accumulated steps so far
  const G4ParticleDefinition* pd; // Particle type (should match all steps)
  G4double Edep;		  // Sum of GetTotalEnergyDeposit()
  G4double Eniel;		  // Sum of GetNonIonizingEnergyDeposit()
  G4ThreeVector start;		  // PreStep position of first hit
  G4ThreeVector end;		  // PostStep position of last hit
};

inline 
std::ostream& operator<<(std::ostream& os, const G4CMPStepAccumulator& accum) {
  accum.Print(os);
  return os;
}

#endif	/* G4CMPStepAccumulator_hh */
