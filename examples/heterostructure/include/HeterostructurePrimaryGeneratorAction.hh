/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructurePrimaryGeneratorAction.hh
/// \brief Definition of the HeterostructurePrimaryGeneratorAction class

// 20260815  Selby Q. Dang (Stanford/SLAC)

#ifndef HeterostructurePrimaryGeneratorAction_h
#define HeterostructurePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"


class G4ParticleGun;
class G4Event;

class HeterostructurePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  HeterostructurePrimaryGeneratorAction();    
  virtual ~HeterostructurePrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                fParticleGun;

};


#endif


