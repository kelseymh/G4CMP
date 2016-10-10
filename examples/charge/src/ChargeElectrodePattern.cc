/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 1e5105a78431cb2647f7bdb285e2c107fa8b1d6a $
//
/// \file  examples/charge/src/ChargeElectrodePattern.cc
/// \brief Simple class demonstrating circumferential electrodes
//
// 20160904  M. Kelsey

#include "ChargeElectrodePattern.hh"
#include "G4AffineTransform.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <cmath>


// Constructor and destructor

ChargeElectrodePattern::ChargeElectrodePattern()
  : G4CMPVElectrodePattern() {;}


// Check coordinates of step, concentric electrodes every 1 cm
G4bool ChargeElectrodePattern::IsNearElectrode(const G4Step& aStep) const {
  if (verboseLevel) 
    G4cout << "ChargeElectrodePattern::IsNearElectrode" << G4endl;

  // Need coordinate transform to get _local_ coordinates on crystal
  // NOTE:  Using pre-step because at boundary PostStep is next volume
  const G4VTouchable* touch = aStep.GetPreStepPoint()->GetTouchable();
  G4AffineTransform toGlobal(touch->GetRotation(), touch->GetTranslation());
  G4AffineTransform toLocal = toGlobal.Inverse();

  G4ThreeVector pos = aStep.GetPreStepPoint()->GetPosition();
  toLocal.ApplyPointTransform(pos);

  // Electrodes are .1 mm rings spaced 2 mm apart on top and bottom faces
  G4double r = pos.rho()/(2.*mm);

  G4bool isnear = ((r - std::floor(r)) < 0.05);

#ifdef G4CMP_DEBUG
  G4cout << " " << aStep.GetTrack()->GetParticleDefinition()->GetParticleName()
	 << " electrode r " << r << " " << (isnear?"yes":"no") << G4endl;
#endif

  return isnear;
}

// Simple absorption to make a hit
void ChargeElectrodePattern::
AbsorbAtElectrode(const G4Track& aTrack, const G4Step& aStep,
		  G4ParticleChange& aParticleChange) const {
  if (verboseLevel) {
    G4cout << "ChargeElectrodePattern::AbsorbAtElectrode "
	   << " rho " << aStep.GetPreStepPoint()->GetPosition().rho()
	   << G4endl;
  }

  //*** FIXME: Need ProcessUtils::GetKineticEnergy() ***
  G4double ekin = aTrack.GetKineticEnergy();
  aParticleChange.ProposeNonIonizingEnergyDeposit(ekin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
}
