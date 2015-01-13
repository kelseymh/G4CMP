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
// 20140312  Introduce multiple inheritance from G4CMPProcessUtils
// 20140325  Move time-step calculation here from TimeStepper and LukeScat
// 20140331  Add required subtype code to constructor
// 20140902  Add new kinematics function which takes energy as input
// 20141231  Add function to enforce minimum step length (fraction of L0)
// 20150112  Rename SetNewKinematics to FillParticleChange for clarity

#ifndef G4CMPVDriftProcess_h
#define G4CMPVDriftProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4CMPProcessSubType.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"


class G4CMPVDriftProcess : public G4VDiscreteProcess, public G4CMPProcessUtils {
public:
  G4CMPVDriftProcess(const G4String& processName, G4CMPProcessSubType stype);
  virtual ~G4CMPVDriftProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

  // Overload base version to set a minimum step size, avoiding "stuck" tracks
  virtual G4double
  PostStepGetPhysicalInteractionLength(const G4Track& track,
				       G4double previousStepSize,
				       G4ForceCondition* condition);

  // NOTE:  These functions must call back to base class implementations!
  virtual void LoadDataForTrack(const G4Track* track);

  // Initialize current valley for currently active track(s)
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

protected:
  // Convenient parameters for computing carrier propagation
  G4double velLong;		// Sound velocity in cystal
  G4double mc_e;		// Effective electron mass
  G4double l0_e;		// Electron scattering length
  G4double ksound_e;		// Electron momentum at sound speed
  G4double mc_h;		// Effective hole mass
  G4double l0_h;		// Hole scattering length
  G4double ksound_h;		// Hole momentum at sound speed

  // Compute characteristic time step for charge carrier
  // Parameters are "Mach number" (ratio with sound speed) and scattering length
  G4double ChargeCarrierTimeStep(G4double mach, G4double l0) const;

  // Minimum Time Step in electric field. Coeff is determined empirically
  G4double ComputeMinTimeStep(const G4Track& track);
  G4double MaxMachStep(G4double Emag, G4double coeff, G4double l0) const;

  // Fill ParticleChange energy and mass for charge carrier of given momentum
  void FillParticleChange(G4int ivalley, const G4ThreeVector& p);

  // Fill ParticleChange energy and mass for charge carrier of given energy
  void FillParticleChange(G4int ivalley, G4double Ekin,
			  const G4ThreeVector& pdir);

private:
  // hide assignment operators as private 
  G4CMPVDriftProcess(G4CMPVDriftProcess&);
  G4CMPVDriftProcess& operator=(const G4CMPVDriftProcess& right);
};

#endif	/* G4CMPVDriftProcess_h */
