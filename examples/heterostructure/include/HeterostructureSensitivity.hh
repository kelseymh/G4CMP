/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructureSensitivity.hh
/// \brief Definition of the HeterostructureSensitivity class

// 20260815  Selby Q. Dang (Stanford/SLAC)

#ifndef HeterostructureSensitivity_h
#define HeterostructureSensitivity_h 1

#include "G4CMPElectrodeHit.hh"
#include "G4CMPElectrodeSensitivity.hh"

class HeterostructureSensitivity final : public G4CMPElectrodeSensitivity {
public:
  HeterostructureSensitivity(G4String);
  virtual ~HeterostructureSensitivity();

  // No copies
  HeterostructureSensitivity(const HeterostructureSensitivity&) = delete;
  HeterostructureSensitivity& operator=(const HeterostructureSensitivity&) = delete;
  HeterostructureSensitivity(HeterostructureSensitivity&&) = delete;
  HeterostructureSensitivity& operator=(HeterostructureSensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream output;
  G4String fileName;
};

#endif  /* HeterostructureSensitivity_h */