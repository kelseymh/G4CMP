/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef PhononSensitivity_h
#define PhononSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class PhononSensitivity final : public G4CMPElectrodeSensitivity {
public:
  PhononSensitivity(G4String name);
  virtual ~PhononSensitivity();
  // No copies
  PhononSensitivity(const PhononSensitivity&) = delete;
  PhononSensitivity& operator=(const PhononSensitivity&) = delete;
  // Move OK
  PhononSensitivity(PhononSensitivity&&);
  PhononSensitivity& operator=(PhononSensitivity&&);

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream output;
  G4String fileName;
};

#endif
