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
/// \file library/include/G4CMPDriftBoundaryProcess.hh
/// \brief Definition of the G4CMPDriftBoundaryProcess base class
//
// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement

#ifndef G4CMPDriftBoundaryProcess_h
#define G4CMPDriftBoundaryProcess_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPDriftBoundaryProcess : public G4CMPVDriftProcess {
public:
  G4CMPDriftBoundaryProcess();
  virtual ~G4CMPDriftBoundaryProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  G4CMPDriftBoundaryProcess(G4CMPDriftBoundaryProcess&);
  G4CMPDriftBoundaryProcess& operator=(const G4CMPDriftBoundaryProcess& right);

  G4double kCarTolerance;
};

#endif
