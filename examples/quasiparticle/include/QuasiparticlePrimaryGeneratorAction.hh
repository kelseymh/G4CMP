/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/quasiparticle/include/QuasiparticlePrimaryGeneratorAction.hh
/// \brief Definition of the QuasiparticlePrimaryGeneratorAction class
//
// $Id: ecbf57649dfaeb88e0fac25491bf8fb68c9308ec $
//

#ifndef QuasiparticlePrimaryGeneratorAction_h
#define QuasiparticlePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"


class G4ParticleGun;
class G4Event;

class QuasiparticlePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  QuasiparticlePrimaryGeneratorAction();    
  virtual ~QuasiparticlePrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4GeneralParticleSource*                fParticleGun;

};


#endif


