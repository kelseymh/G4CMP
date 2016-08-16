/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef PhononSensitivity_h
#define PhononSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class G4HCofThisEvent;

using std::ofstream;

class PhononSensitivity : public G4CMPElectrodeSensitivity
{
public:
  PhononSensitivity(G4String name);
  virtual ~PhononSensitivity();
  //virtual void Initialize(G4HCofThisEvent*);
  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);

private:
  ofstream output;
  G4String fileName;
};

#endif
