#ifndef ChargeElectrodeSensitivity_h
#define ChargeElectrodeSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPElectrodeHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class ofstream;

using std::ofstream;

class ChargeElectrodeSensitivity
    : public G4CMPElectrodeSensitivity
{
public:
  ChargeElectrodeSensitivity(G4String);
  virtual ~ChargeElectrodeSensitivity();
  
  virtual void EndOfEvent(G4HCofThisEvent*);
  
private:
  ofstream output;
};

#endif
