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

#ifndef G4CMPeLukeScattering_h
#define G4CMPeLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4VProcess;


class G4CMPeLukeScattering : public G4CMPVDriftProcess {
public:
  G4double MakeTheta(G4double& k,G4double& ks);
  G4double MakePhi(G4double& k, G4double& ks, G4double& theta);

  G4CMPeLukeScattering(G4VProcess* stepper);

  virtual ~G4CMPeLukeScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  virtual bool IsApplicable(const G4ParticleDefinition&);

protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

private:
  //hide assignment operator as private
  G4CMPeLukeScattering(G4CMPeLukeScattering&);
  G4CMPeLukeScattering& operator=(const G4CMPeLukeScattering& right);

  G4VProcess* stepLimiter;

  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;
  G4RotationMatrix mInv; 	// Inverse mass tensor
  G4ThreeVector T;       	// HV Transformation matrix diagonal
};

#endif
