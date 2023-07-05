/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononPrimaryGeneratorAction.hh
/// \brief Definition of the PhononPrimaryGeneratorAction class
//
// $Id: ecbf57649dfaeb88e0fac25491bf8fb68c9308ec $
//

#ifndef PhononPrimaryGeneratorAction_h
#define PhononPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"


class G4ParticleGun;
class G4Event;

class PhononPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PhononPrimaryGeneratorAction();    
  virtual ~PhononPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                fParticleGun;

};


#endif


