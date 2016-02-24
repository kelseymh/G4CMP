/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$

#include "G4CMPhDriftBoundaryProcess.hh"
#include "G4CMPDriftHole.hh"
#include "G4LatticeManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4CMPhDriftBoundaryProcess::G4CMPhDriftBoundaryProcess()
  : G4CMPVDriftBoundaryProcess("hDrift", G4CMPDriftHole::Definition()) {;}

G4CMPhDriftBoundaryProcess::~G4CMPhDriftBoundaryProcess() {}

// Apply kinematic absoprtion (wave-vector at surface)

G4bool G4CMPhDriftBoundaryProcess::AbsorbTrack(const G4Step& aStep) {
  return (G4CMPVDriftBoundaryProcess::AbsorbTrack(aStep) ||
    GetLocalWaveVector(*(aStep.GetTrack()))*surfNorm > absMinKHole);
}
