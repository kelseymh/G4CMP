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
/// \file library/include/G4CMPeDriftBoundaryProcess.hh
/// \brief Definition of the G4CMPeDriftBoundaryProcess class
//
// $Id$
//
// 20150601  M. Kelsey -- Follow encapsulation of boundary process actions

#ifndef G4CMPeDriftBoundaryProcess_h
#define G4CMPeDriftBoundaryProcess_h 1

#include "globals.hh"
#include "G4CMPVDriftBoundaryProcess.hh"

class G4CMPeDriftBoundaryProcess : public G4CMPVDriftBoundaryProcess {
public:
  G4CMPeDriftBoundaryProcess();
  virtual ~G4CMPeDriftBoundaryProcess();

protected:
  virtual G4ThreeVector GetLocalWaveVector(const G4Track& aTrack) const;

  // Apply kinematic absoprtion (wave-vector at surface)
  virtual G4bool AbsorbTrack(const G4Step&);

  // Apply reflection to velocity vector, not momentum
  virtual G4bool ReflectTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoReflection(const G4Step& aStep);

private:
  G4CMPeDriftBoundaryProcess(G4CMPeDriftBoundaryProcess&);
  G4CMPeDriftBoundaryProcess& operator=(const G4CMPeDriftBoundaryProcess& right);
};

#endif
