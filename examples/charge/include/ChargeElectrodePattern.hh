/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 1e5105a78431cb2647f7bdb285e2c107fa8b1d6a $
//
/// \file  examples/charge/include/ChargeElectrodePattern.hh
/// \brief Simple class demonstrating circumferential electrodes
//
// 20160904  M. Kelsey
// 20200601  G4CMP-207:  Implement clone function for thread-local copying

#ifndef ChargeElectrodePattern_h
#define ChargeElectrodePattern_h 1

#include "G4CMPVElectrodePattern.hh"


class ChargeElectrodePattern : public G4CMPVElectrodePattern {
public:
  ChargeElectrodePattern();
  virtual ~ChargeElectrodePattern() {;}

  // Use default copy/move operators
  ChargeElectrodePattern(const ChargeElectrodePattern&) = default;
  ChargeElectrodePattern(ChargeElectrodePattern&&) = default;
  ChargeElectrodePattern& operator=(const ChargeElectrodePattern&) = default;
  ChargeElectrodePattern& operator=(ChargeElectrodePattern&&) = default;

  virtual G4CMPVElectrodePattern* Clone() const {
    return new ChargeElectrodePattern(*this);
  }

  virtual G4bool IsNearElectrode(const G4Step& aStep) const;
  virtual void AbsorbAtElectrode(const G4Track& aTrack, const G4Step& aStep,
				 G4ParticleChange& aParticleChange) const;
};

#endif	/* ChargeElectrodePattern_h */
