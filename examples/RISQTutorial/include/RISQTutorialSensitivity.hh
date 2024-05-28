/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialSensitivity_h
#define RISQTutorialSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class RISQTutorialSensitivity final : public G4CMPElectrodeSensitivity {
public:
  RISQTutorialSensitivity(G4String name);
  virtual ~RISQTutorialSensitivity();
  // No copies
  RISQTutorialSensitivity(const RISQTutorialSensitivity&) = delete;
  RISQTutorialSensitivity& operator=(const RISQTutorialSensitivity&) = delete;
  /* Move is disabled for now because old versions of GCC can't move ofstream
  // Move OK
  RISQTutorialSensitivity(RISQTutorialSensitivity&&);
  RISQTutorialSensitivity& operator=(RISQTutorialSensitivity&&);
  */
  RISQTutorialSensitivity(RISQTutorialSensitivity&&) = delete;
  RISQTutorialSensitivity& operator=(RISQTutorialSensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetHitOutputFile(const G4String& fn);
  void SetPrimaryOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream primaryOutput;
  std::ofstream hitOutput;
  G4String primaryFileName;
  G4String hitFileName;
  
};

#endif
