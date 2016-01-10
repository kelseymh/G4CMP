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
/// \file library/src/G4CMPeLukeScattering.cc
/// \brief Implementation of the G4CMPeLukeScattering class
//
// $Id$
//
// 20140325  Move time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140415  Add run-time flag to select valley vs. H-V kinematics
// 20140430  Compute kinematics using mass tensor; prepare to create phonons
// 20140509  Remove valley vs. H-V flag; add run-time envvar to bias phonons
// 20140521  Remove momentum-check loop; energy conservation is enforced
// 20140903  Get Etrack using valley kinematics, _not_ track or stepPoint
// 201411??  R.Agnese -- Merge functionality from TimeStepper here
// 20141231  Rename "minimum step" function to ComputeMinTimeStep
// 20150106  Move envvar to G4CMPConfigManager, unify overlapping code w/hLuke
// 20150111  Migrate most physics to new base class G4CMPVLukeScattering,
//	     follow renaming of SetNewKinematics to FillParticleChange

#include "G4CMPeLukeScattering.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
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

G4CMPeLukeScattering::G4CMPeLukeScattering(G4VProcess* stepper)
  : G4CMPVLukeScattering("eLuke", G4CMPDriftElectron::Definition(), stepper) {;}

G4CMPeLukeScattering::~G4CMPeLukeScattering() {;}


// Physics

G4ThreeVector G4CMPeLukeScattering::GetLocalWaveVector(const G4Track& aTrack) const {
  return theLattice->MapV_elToK_HV(GetValleyIndex(aTrack),
                                   GetLocalVelocityVector(aTrack));
}

// Convert local wave-vector to global using HV transform

void G4CMPeLukeScattering::MakeGlobalPhonon(G4ThreeVector& kphonon) const {
  kphonon = theLattice->MapK_HVtoK_valley(GetValleyIndex(GetCurrentTrack()),kphonon);
  RotateToGlobalDirection(kphonon);
}

// Convert local wave-vector to global momentum using HV transform

void G4CMPeLukeScattering::MakeGlobalRecoil(G4ThreeVector& krecoil) const {
  krecoil = theLattice->MapK_HVtoP(GetValleyIndex(GetCurrentTrack()), krecoil);
  RotateToGlobalDirection(krecoil);
}
