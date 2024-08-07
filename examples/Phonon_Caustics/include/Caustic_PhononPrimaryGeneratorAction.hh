/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/


//20240110 Israel Hernandez -- Illinois Institute of Technology, Quantum Science Center and Fermilab


#ifndef Caustic_PhononPrimaryGeneratorAction_h
#define Caustic_PhononPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"


class G4ParticleGun;
class G4Event;

class Caustic_PhononPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Caustic_PhononPrimaryGeneratorAction();
  virtual ~Caustic_PhononPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
  
    G4GeneralParticleSource *fParticleGun;

};


#endif
