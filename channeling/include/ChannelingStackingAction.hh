/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file hadronic/Hadr01/include/ChannelingStackingAction.hh
/// \brief Definition of the ChannelingStackingAction class
//
// $Id: 1508c59183c98bb83f29fbdf81038b1ef036608b $
//
/////////////////////////////////////////////////////////////////////////
//
// ChannelingStackingAction
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef ChannelingStackingAction_h
#define ChannelingStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

class ChannelingStackingMessenger;
class G4Track;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChannelingStackingAction : public G4UserStackingAction
{
public:

  ChannelingStackingAction();
  virtual ~ChannelingStackingAction();
   
  void SetKillStatus(G4bool value);
  void SetKill(const G4String& name);
     
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    
private:

  ChannelingStackingMessenger*  fStackMessenger;
  G4bool              fKillSecondary;

  const G4ParticleDefinition* fParticle;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

