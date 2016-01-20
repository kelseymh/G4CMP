#ifndef ChargeElectrodeSensitivity_h
#define ChargeElectrodeSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class G4HCofThisEvent;
class ChargeFETDigitizerModule;

using std::ofstream;

class ChargeElectrodeSensitivity
    : public G4CMPElectrodeSensitivity
{
public:
  ChargeElectrodeSensitivity(G4String name);
  virtual ~ChargeElectrodeSensitivity();
  //virtual void Initialize(G4HCofThisEvent*);
  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetOutputFile(const G4String& fn);
  
private:
  ofstream output;
  G4String fileName;
  ChargeFETDigitizerModule* FET;
};

#endif
