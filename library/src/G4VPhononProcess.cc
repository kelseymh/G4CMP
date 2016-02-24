/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4VPhononProcess.cc
/// \brief Implementation of the G4VPhononProcess base class
//
// $Id$
//
// 20131111  Add verbosity to report creating secondaries
// 20140312  Move utility functions to new G4CMPProcessUtils
// 20140331  Add required process subtype code

#include "G4VPhononProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4ProcessType.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

namespace {
  const G4ThreeVector nullVec(0.,0.,0.);	// For convenience below
}

// Constructor and destructor

G4VPhononProcess::G4VPhononProcess(const G4String& processName,
				   G4CMPProcessSubType stype)
  : G4VDiscreteProcess(processName, fPhonon), G4CMPProcessUtils() {
  SetProcessSubType(stype);
}

G4VPhononProcess::~G4VPhononProcess() {;}


// Only applies to the known phonon polarization states

G4bool G4VPhononProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4PhononLong::Definition() ||
	  &aPD==G4PhononTransFast::Definition() ||
	  &aPD==G4PhononTransSlow::Definition() );
}


// Initialize wave vectors for currently active track(s)

void G4VPhononProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions
  LoadDataForTrack(track);
}

void G4VPhononProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  ReleaseTrack();
}
