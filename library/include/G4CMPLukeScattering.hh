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
/// \file library/include/G4CMPLukeScattering.hh
/// \brief Definition of the G4CMPLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20160110  Remerge the electron and hole subclasses into one class

#ifndef G4CMPLukeScattering_h
#define G4CMPLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4ThreeVector.hh"
#include <iostream>

class G4CMPTrackInformation;
class G4VProcess;
class G4ParticleDefinition;
class G4Track;

class G4CMPLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPLukeScattering(G4VProcess* stepper=0);
  virtual ~G4CMPLukeScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4double CalculateKSound(const G4CMPTrackInformation* trackInfo);
  // hide assignment operator as private
  G4CMPLukeScattering(G4CMPLukeScattering&);
  G4CMPLukeScattering& operator=(const G4CMPLukeScattering& right);

  G4VProcess* stepLimiter;
#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4CMPLukeScattering */
