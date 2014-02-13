#ifndef G4CMPElectrodeSensitivity_h
#define G4CMPElectrodeSensitivity_h 1

#include "G4VSensitiveDetector.hh"
#include "G4CMPElectrodeHit.hh"

#include <iostream>
#include <fstream>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

using namespace std;


class G4CMPElectrodeSensitivity : public G4VSensitiveDetector {
public:
  G4CMPElectrodeSensitivity(G4String);
  virtual ~G4CMPElectrodeSensitivity();
  
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  
  G4CMPElectrodeHitsCollection* getHitsCollection();
  static G4CMPElectrodeHitsCollection* hitsCollection;
  
private:
  ofstream writer; //writing hit posn to file. Temporary fix.
  ofstream writer2; //writing timing information to file. Temporary fix.
  
  G4int HCID;
};

#endif
