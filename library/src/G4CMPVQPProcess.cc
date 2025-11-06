/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVQPProcess.cc
/// \brief Implementation of the G4CMPVQPProcess base class
//
// $Id$
//
// 20131111  Add verbosity to report creating secondaries
// 20140312  Move utility functions to new G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20170601  Inherit from new G4CMPVProcess, only need IsApplicable

#include "G4CMPVQPProcess.hh"
#include "G4CMPUtils.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4RandomDirection.hh"

// Only applies to QPs

G4bool G4CMPVQPProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsQP(&aPD);
}

//This standardizes the directionality of the final state momenta for all QP
//particles, assuming an XY-oriented film. This should be run at the end
//of all postStepDoIts for all purely discrete QP processes where the QPs
//still exist
G4bool G4CMPVQPProcess::RandomizeFinalStateMomentumDirectionInXY() {
  G4ThreeVector returnDir = G4RandomDirection();
  returnDir.setZ(0);
  aParticleChange.ProposeMomentumDirection(returnDir.unit());
  return true;
}
