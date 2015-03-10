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
/// \file library/include/G4CMPIonisationWrapper.hh
/// \brief Definition of the G4CMPIonisationWrapper process class.  This
///	class will be used to extend the existing Geant4 ionization
///	(and possibly other) processes to generate phonons and charge
///	carriers as secondaries.
//
// $Id$
//
// 20150306  Michael Kelsey

#ifndef G4CMPIonisationWrapper_hh
#define G4CMPIonisationWrapper_hh 1

#include "G4WrapperProcess.hh"
#include "G4CMPProcessUtils.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <utility>

class G4Track;
class G4Step;
class G4VParticleChange;


class G4CMPIonisationWrapper : public G4WrapperProcess,
			       public G4CMPProcessUtils {
public:
  G4CMPIonisationWrapper(G4VProcess* ionizProc);
  virtual ~G4CMPIonisationWrapper();

  // Call through to wrapped process to get energy loss
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
					  const G4Step& stepData);

  virtual G4VParticleChange* AlongStepDoIt(const G4Track& track,
					   const G4Step& stepData);

protected:
  void AddPhonons(G4VParticleChange* theChange, const G4Step& stepData);

  void AddChargeCarriers(G4VParticleChange* theChange, const G4Step& stepData);

  void GenerateEnergyPositions(const G4Step& stepData, G4double Etotal,
			       G4double Eunit, G4double sigmaE);

private:
  G4double energyPerPhonon;	// Energy partition and fluctuations
  G4double sigmaEPerPhonon;

  G4double energyPerChargePair;
  G4double sigmaEPerChargePair;

  typedef std::pair<G4double, G4double> EPosPair;
  std::vector<EPosPair> energyPosList;		// Buffer for secondaries

  // No copying allowed
  G4CMPIonisationWrapper(const G4CMPIonisationWrapper& right);
  G4CMPIonisationWrapper& operator=(const G4CMPIonisationWrapper& right);
};

#endif	/* G4CMPIonisationWrapper_hh */
