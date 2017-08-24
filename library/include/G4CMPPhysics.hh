/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// Create particles and physics processes for phonons and charge carriers
// Usage:  [physics-list]->AddPhysics(new G4CMPPhysics);
//
// 20150309  M. Kelsey -- Add function to find and wrap *Ionisation processes

#ifndef G4CMPPhysics_hh
#define G4CMPPhysics_hh 1

#include "G4VPhysicsConstructor.hh"


class G4CMPPhysics : public G4VPhysicsConstructor {
public:
  G4CMPPhysics(const G4String& name="G4CMPPhysics");
  virtual ~G4CMPPhysics() {;}
  
public:
  virtual void ConstructParticle();	// Creates phonons and "drifters"
  virtual void ConstructProcess();	// Adds processes to physics list

protected:
  void AddSecondaryProduction();	// All charged particles make e/h, phn

private:
  G4CMPPhysics(const G4CMPPhysics& rhs);		// Copying is forbidden
  G4CMPPhysics& operator=(const G4CMPPhysics& rhs);
};

#endif	/* G4CMPPhysics_hh */
