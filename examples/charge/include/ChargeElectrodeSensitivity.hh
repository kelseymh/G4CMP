/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// 20170830  Remove FET simulation

#ifndef ChargeElectrodeSensitivity_h
#define ChargeElectrodeSensitivity_h 1

#include "G4CMPElectrodeHit.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include <memory>


class ChargeElectrodeSensitivity final : public G4CMPElectrodeSensitivity {
public:
  ChargeElectrodeSensitivity(G4String);
  virtual ~ChargeElectrodeSensitivity();

  // No copies
  ChargeElectrodeSensitivity(const ChargeElectrodeSensitivity&) = delete;
  ChargeElectrodeSensitivity& operator=(const ChargeElectrodeSensitivity&) = delete;
  ChargeElectrodeSensitivity(ChargeElectrodeSensitivity&&) = delete;
  ChargeElectrodeSensitivity& operator=(ChargeElectrodeSensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream output;
  G4String fileName;
};

#endif
