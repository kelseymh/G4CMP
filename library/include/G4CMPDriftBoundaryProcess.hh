/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPDriftBoundaryProcess.hh
/// \brief Definition of the G4CMPDriftBoundaryProcess base class
//
// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20150212  Remove file IO. Use sensitive detectors instead
// 20150304  Change to generic G4CMPVDrifBoundaryProcess and
//           utilize specific G4CMPDrift{Electron,Hole}BoundaryProcess
// 20150420  Replace MFP with GPIL to suppress unnecessary verbosity.
// 20150529  Add DoReflection() function, so electron can overload
// 20150603  Add parameter to count number of reflections by track
// 20160215  Reunify boundary processes at least until we remove CDMS-specific
//           stuff.

#ifndef G4CMPDriftBoundaryProcess_h
#define G4CMPDriftBoundaryProcess_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPDriftBoundaryProcess : public G4CMPVDriftProcess {
public:
  G4CMPDriftBoundaryProcess(const G4String& name="Drift");
  virtual ~G4CMPDriftBoundaryProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                   G4double previousStepSize,
                                                   G4ForceCondition* condition);

  // NOTE:  These functions must call back to base class implementations!
  virtual void LoadDataForTrack(const G4Track* track);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoAbsorption(const G4Step& aStep);

  virtual G4bool ReflectTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoReflection(const G4Step& aStep);

private:
  G4CMPDriftBoundaryProcess(G4CMPDriftBoundaryProcess&);
  G4CMPDriftBoundaryProcess& operator=(const G4CMPDriftBoundaryProcess& right);

protected:
  G4double kCarTolerance;
  G4String shortName;
  const G4ParticleDefinition* theCarrier;

  G4double reflProb;
  G4double absProb;
  G4double absMinK;
  G4double maxRefl;
  G4ThreeVector surfNorm;	// Surface normal (temporary buffer)

  G4int numberOfReflections;	// Counter to prevent runaway tracks
};

#endif
