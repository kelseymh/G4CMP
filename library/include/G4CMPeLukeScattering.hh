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
/// \file library/include/G4eLukeScattering.hh
/// \brief Definition of the G4eLukeScattering class
//
// $Id$
//
// 20140404  Drop unnecessary data members, using functions in G4LatticePhysical
// 20140407  Move angle generation to base class
// 20150111  Add debugging output file with phonon energies
// 20150111  Move majority of interface to new base class

#ifndef G4CMPeLukeScattering_h
#define G4CMPeLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVLukeScattering.hh"

class G4VProcess;


class G4CMPeLukeScattering : public G4CMPVLukeScattering {
public:
  G4CMPeLukeScattering(G4VProcess* stepper);
  virtual ~G4CMPeLukeScattering();

  //*****  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4ThreeVector GetLocalWaveVector(const G4Track& aTrack) const;

  // Convert local wave-vector to global using HV transform
  virtual void MakeGlobalPhonon(G4ThreeVector& kphonon) const;

  // Convert local wave-vector to global momentum using HV transform
  virtual void MakeGlobalRecoil(G4ThreeVector& krecoil) const;

private:
  //hide assignment operator as private
  G4CMPeLukeScattering(G4CMPeLukeScattering&);
  G4CMPeLukeScattering& operator=(const G4CMPeLukeScattering& right);
};

#endif
