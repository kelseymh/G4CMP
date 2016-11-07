/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChargeElectrodeSensitivity_h
#define ChargeElectrodeSensitivity_h 1

#include "G4CMPElectrodeHit.hh"
#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class ChargeFETDigitizerModule;

class ChargeElectrodeSensitivity : public G4VSensitiveDetector {
public:
  ChargeElectrodeSensitivity(G4String name);
  virtual ~ChargeElectrodeSensitivity() override;
  // No copies
  ChargeElectrodeSensitivity(const ChargeElectrodeSensitivity&) = delete;
  ChargeElectrodeSensitivity& operator=(const ChargeElectrodeSensitivity&) = delete;
  // Move OK
  ChargeElectrodeSensitivity(ChargeElectrodeSensitivity&&);
  ChargeElectrodeSensitivity& operator=(ChargeElectrodeSensitivity&&);

  virtual void Initialize(G4HCofThisEvent*) override;
  virtual void EndOfEvent(G4HCofThisEvent*) override;

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*) override;

private:
  G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;
  std::unique_ptr<ChargeFETDigitizerModule> FET;
  G4CMPElectrodeHitsCollection* hitsCollection;
  std::ofstream output;
  G4String fileName;
};

#endif
