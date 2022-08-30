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
// 20170601  Inherit from new G4CMPVProcess, only need IsApplicable

#include "G4VPhononProcess.hh"
#include "G4CMPUtils.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"


// Only applies to the known phonon polarization states

G4bool G4VPhononProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsPhonon(&aPD);
}
