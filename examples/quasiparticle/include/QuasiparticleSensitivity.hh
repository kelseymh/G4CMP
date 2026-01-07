/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleSensitivity_h
#define QuasiparticleSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class QuasiparticleSensitivity final : public G4CMPElectrodeSensitivity {
public:
  QuasiparticleSensitivity(G4String name);
  virtual ~QuasiparticleSensitivity();
  // No copies
  QuasiparticleSensitivity(const QuasiparticleSensitivity&) = delete;
  QuasiparticleSensitivity& operator=(const QuasiparticleSensitivity&) = delete;
  /* Move is disabled for now because old versions of GCC can't move ofstream
  // Move OK
  QuasiparticleSensitivity(QuasiparticleSensitivity&&);
  QuasiparticleSensitivity& operator=(QuasiparticleSensitivity&&);
  */
  QuasiparticleSensitivity(QuasiparticleSensitivity&&) = delete;
  QuasiparticleSensitivity& operator=(QuasiparticleSensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream output;
  G4String fileName;
};

#endif
