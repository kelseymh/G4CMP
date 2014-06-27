//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr01/include/ChannelingStackingAction.hh
/// \brief Definition of the ChannelingStackingAction class
//
// $Id$
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

