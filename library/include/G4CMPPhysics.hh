// $Id$
//
// Create particles and physics processes for phonons and charge carriers
// Usage:  [physics-list]->AddPhysics(new G4CMPPhysics);

#ifndef G4CMPPhysics_hh
#define G4CMPPhysics_hh 1

#include "G4VPhysicsConstructor.hh"


class G4CMPPhysics : public G4VPhysicsConstructor {
public:
  G4CMPPhysics(const G4String& name="G4CMPPhysics")
    : G4VPhysicsConstructor(name) {;}

  virtual ~G4CMPPhysics() {;}

public:
  virtual void ConstructParticle();	// Creates phonons and "drifters"
  virtual void ConstructProcess();	// Adds processes to physics list

private:
  G4CMPPhysics(const G4CMPPhysics& rhs);		// Copying is forbidden
  G4CMPPhysics& operator=(const G4CMPPhysics& rhs);
};

#endif	/* G4CMPPhysics_hh */
