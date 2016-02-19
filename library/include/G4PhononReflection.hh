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
/// \file library/include/G4PhononReflection.hh
/// \brief Definition of the G4PhononReflection class
//
// $Id$
//
#ifndef G4PhononReflection_h
#define G4PhononReflection_h 1

#include "G4VPhononProcess.hh"

class G4CMPTrackInformation;
class G4CMPSurfaceProperty;

class G4PhononReflection : public G4VPhononProcess {
public:
  G4PhononReflection(const G4String& processName ="phononReflection" );
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep);
  
protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double prevStepLength,
                                   G4ForceCondition* condition);
  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoAbsorption(const G4Step& aStep);

  virtual G4bool ReflectTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoReflection(const G4Step& aStep);

  G4bool ReflectionIsGood(G4int polarization);

  G4double kCarTolerance;
  G4double absProb;
  G4double reflProb;
  G4double specProb;
  G4double absMinK;
  G4double maxRefl;
  G4CMPTrackInformation* trackInfo;
  G4ThreeVector waveVector;
  G4ThreeVector surfNorm;

private:
  // hide assignment operator as private
  G4PhononReflection(G4PhononReflection&);
  G4PhononReflection& operator=(const G4PhononReflection& right);
};

#endif










