/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef ChargePrimaryGeneratorAction_h
#define ChargePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;

class ChargePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ChargePrimaryGeneratorAction();    
  virtual ~ChargePrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                particleGun;

};


#endif


