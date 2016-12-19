/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPElectrodeSensitivity_hh
#define G4CMPElectrodeSensitivity_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4CMPElectrodeHit.hh"

class G4HCofThisEvent;

class G4CMPElectrodeSensitivity : public G4VSensitiveDetector {
public:
  explicit G4CMPElectrodeSensitivity(G4String);
  // No copies
  G4CMPElectrodeSensitivity(const G4CMPElectrodeSensitivity&) = delete;
  G4CMPElectrodeSensitivity& operator=(const G4CMPElectrodeSensitivity&) = delete;
  // Move OK
  G4CMPElectrodeSensitivity(G4CMPElectrodeSensitivity&&);
  G4CMPElectrodeSensitivity& operator=(G4CMPElectrodeSensitivity&&);

  virtual void Initialize(G4HCofThisEvent*) override;
  
protected:
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

  G4CMPElectrodeHitsCollection* hitsCollection;
};

#endif
