/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/RISQTutorialPrimaryGeneratorAction.hh
/// \brief Definition of the RISQTutorialPrimaryGeneratorAction class
//
// $Id: ecbf57649dfaeb88e0fac25491bf8fb68c9308ec $
//

#ifndef RISQTutorialPrimaryGeneratorAction_h
#define RISQTutorialPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"


class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

class RISQTutorialPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  RISQTutorialPrimaryGeneratorAction();    
  virtual ~RISQTutorialPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4GeneralParticleSource*                fParticleGun;

};


#endif


