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
/// \file library/include/G4CMPVDriftBoundaryProcess.hh
/// \brief Definition of the G4CMPVDriftBoundaryProcess base class
//
// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils
// 20140331  Inherit from G4CMPVDriftProcess to get subtype enforcement
// 20150212  Remove file IO. Use sensitive detectors instead
// 20150304  Change to generic G4CMPVDrifBoundaryProcess and
//           utilize specific G4CMPDrift{Electron,Hole}BoundaryProcess
// 20150420  Replace MFP with GPIL to suppress unnecessary verbosity.
// 20150529  Add DoReflection() function, so electron can overload

#ifndef G4CMPVDriftBoundaryProcess_h
#define G4CMPVDriftBoundaryProcess_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"


class G4CMPVDriftBoundaryProcess : public G4CMPVDriftProcess {
public:
  G4CMPVDriftBoundaryProcess(const G4String& name="Drift",
               const G4ParticleDefinition* carrier=0);
  virtual ~G4CMPVDriftBoundaryProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& track,
				       G4double previousStepSize,
				       G4ForceCondition* condition);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  // NOTE:  These functions must call back to base class implementations!
  virtual void LoadDataForTrack(const G4Track* track);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
  virtual G4ThreeVector GetWaveVector(const G4Track&) const =0;
  virtual G4double GetKineticEnergy(const G4Track&) const =0;

  // Decide and apply different surface actions; subclasses may override
  virtual G4bool AbsorbTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoAbsorption(const G4Step& aStep);

  virtual G4bool HitElectrode(const G4Step& aStep);
  virtual G4VParticleChange* DoElectrodeHit(const G4Step& aStep);

  virtual G4bool ReflectTrack(const G4Step& aStep);
  virtual G4VParticleChange* DoReflection(const G4Step& aStep);

private:
  G4CMPVDriftBoundaryProcess(G4CMPVDriftBoundaryProcess&);
  G4CMPVDriftBoundaryProcess& operator=(const G4CMPVDriftBoundaryProcess& right);

  G4double kCarTolerance;

protected:
  const G4ParticleDefinition* theCarrier;
  G4String shortName;

  G4double absProb;
  G4double absDeltaV;
  G4double absMinKElec;
  G4double absMinKHole;
  G4double electrodeV;
  G4ThreeVector surfNorm;	// Surface normal (temporary buffer)
};

#endif
