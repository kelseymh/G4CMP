/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


/// \file ValidationSensitivity.hh
/// \brief Definition of the validation example's sensitivity class

#ifndef ValidationSensitivity_h
#define ValidationSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class ValidationSensitivity final : public G4CMPElectrodeSensitivity {
public:
  ValidationSensitivity(G4String name);
  virtual ~ValidationSensitivity();
  // No copies
  ValidationSensitivity(const ValidationSensitivity&) = delete;
  ValidationSensitivity& operator=(const ValidationSensitivity&) = delete;
  /* Move is disabled for now because old versions of GCC can't move ofstream
  // Move OK
  ValidationSensitivity(ValidationSensitivity&&);
  ValidationSensitivity& operator=(ValidationSensitivity&&);
  */
  ValidationSensitivity(ValidationSensitivity&&) = delete;
  ValidationSensitivity& operator=(ValidationSensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream output;
  G4String fileName;
};

#endif
