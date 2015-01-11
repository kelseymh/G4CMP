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
/// \file library/include/G4CMPhLukeScattering.hh
/// \brief Definition of the G4CMPhLukeScattering class
//
// $Id$
//
// 20140407  Move angle generation to base class
// 20150111  Add debugging output file with phonon energies

#ifndef G4CMPhLukeScattering_h
#define G4CMPhLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include <iostream>

class G4VProcess;


class G4CMPhLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPhLukeScattering(G4VProcess* stepper);
  virtual ~G4CMPhLukeScattering();
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
  
  // Override default applicability for phonon processes
  virtual bool IsApplicable(const G4ParticleDefinition&);
  
protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  G4CMPhLukeScattering(G4CMPhLukeScattering&);
  G4CMPhLukeScattering& operator=(const G4CMPhLukeScattering& right);

  G4VProcess* stepLimiter;

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4CMPhLukeScattering */
