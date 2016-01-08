//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file library/src/G4CMPhLukeScattering.cc
/// \brief Implementation of the G4CMPhLukeScattering class
//
// $Id$
//
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140509  Add run-time envvar to bias phonons
// 20141231  Rename "minimum step" function to ComputeMinTimeStep
// 20150111  Move envvar to G4CMPConfigManager, unify overlapping code w/eLuke
// 20150111  Migrate most physics to new base class G4CMPVLukeScattering

#include "G4CMPhLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftHole.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include <iostream>
#include <fstream>


// Constructor and destructor

G4CMPhLukeScattering::G4CMPhLukeScattering(G4VProcess* stepper)
  : G4CMPVLukeScattering("hLuke", G4CMPDriftHole::Definition(), stepper) {;}

G4CMPhLukeScattering::~G4CMPhLukeScattering() {;}


// Physics

G4ThreeVector G4CMPhLukeScattering::GetLocalWaveVector(const G4Track& aTrack) const {
  return GetLocalMomentum(aTrack) / hbarc;
}

// Convert local wave-vector to global

void G4CMPhLukeScattering::MakeGlobalPhonon(G4ThreeVector& kphonon) const {
  RotateToGlobalDirection(kphonon);
}

// Convert local wave-vector to global momentum

void G4CMPhLukeScattering::MakeGlobalRecoil(G4ThreeVector& krecoil) const {
  krecoil *= hbarc;
  RotateToGlobalDirection(krecoil);
}
