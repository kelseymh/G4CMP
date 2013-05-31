#ifndef AlminumElectrodeSensitivity_h
#define AlminumElectrodeSensitivity_h 1

#include "G4VSensitiveDetector.hh"
#include "AlminumElectrodeHit.hh"

#include <iostream>
#include <fstream>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

using namespace std;


class AlminumElectrodeSensitivity : public G4VSensitiveDetector
{

  public:
      AlminumElectrodeSensitivity(G4String);
      virtual ~AlminumElectrodeSensitivity();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);

  AlminumElectrodeHitsCollection* getHitsCollection();
  static AlminumElectrodeHitsCollection* hitsCollection;

  private:
  //AlminumElectrodeHitsCollection * hitsCollection;
  ofstream writer; //writing hit posn to file. Temporary fix.
  ofstream writer2; //writing timing information to file. Temporary fix.

      G4int HCID;
};




#endif

