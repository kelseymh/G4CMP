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
/// \file library/include/G4CMPVDriftProcess.hh
/// \brief Definition of the G4CMPVDriftProcess base class
//
// $Id$
//
#ifndef G4CMPVDriftProcess_h
#define G4CMPVDriftProcess_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"


class G4CMPVDriftProcess : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4CMPVDriftProcess(const G4String& processName);
  virtual ~G4CMPVDriftProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

  // Initialize current valley for currently active track(s)
  // NOTE:  These functions must call back to base class implementations!
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

private:
  // hide assignment operators as private 
  G4CMPVDriftProcess(G4CMPVDriftProcess&);
  G4CMPVDriftProcess& operator=(const G4CMPVDriftProcess& right);
};

#endif	/* G4CMPVDriftProcess_h */
