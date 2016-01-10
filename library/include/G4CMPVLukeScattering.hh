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
/// \file library/include/G4CMPVLukeScattering.hh
/// \brief Definition of the G4CMPVLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes

#ifndef G4CMPVLukeScattering_h
#define G4CMPVLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4ThreeVector.hh"
#include <iostream>

class G4VProcess;
class G4ParticleDefinition;
class G4Track;

class G4CMPVLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPVLukeScattering(const G4String& name="Luke",
		       const G4ParticleDefinition* carrier=0,
		       G4VProcess* stepper=0);
  virtual ~G4CMPVLukeScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  // Override default applicability for phonon processes
  virtual bool IsApplicable(const G4ParticleDefinition&);

  // NOTE:  These functions must call back to base class implementations!
  virtual void LoadDataForTrack(const G4Track* track);
  
protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

  // Charge carrier subclasses must implement these individually
  virtual G4ThreeVector GetLocalWaveVector(const G4Track& aTrack) const = 0;
  virtual G4double GetWaveNumber(const G4Track& aTrack) const;

  // Convert local wave-vector to global (electrons need HV transform)
  virtual void MakeGlobalPhonon(G4ThreeVector& kphonon) const = 0;

  // Convert local wave-vector to global momentum (electrons need HV transform)
  virtual void MakeGlobalRecoil(G4ThreeVector& krecoil) const = 0;

private:
  // hide assignment operator as private
  G4CMPVLukeScattering(G4CMPVLukeScattering&);
  G4CMPVLukeScattering& operator=(const G4CMPVLukeScattering& right);

protected:
  G4String shortName;
  G4VProcess* stepLimiter;
  std::ofstream output;

  // Particle-specific parameters for use by MFP and DoIt
  const G4ParticleDefinition* theCarrier;
  G4double theKsound;
  G4double theL0;
};

#endif	/* G4CMPVLukeScattering */
