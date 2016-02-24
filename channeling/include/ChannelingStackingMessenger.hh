/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file hadronic/Hadr01/include/ChannelingStackingMessenger.hh
/// \brief Definition of the ChannelingStackingMessenger class
//
// $Id: 08e0731fe2ab53d8df28ddd2cbb8dbeed5611cfe $
//
//
/////////////////////////////////////////////////////////////////////////
//
// ChannelingStackingMessenger
//
// Created: 31.05.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

#ifndef ChannelingStackingMessenger_h
#define ChannelingStackingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class ChannelingStackingAction;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChannelingStackingMessenger: public G4UImessenger
{
public:

  ChannelingStackingMessenger(ChannelingStackingAction*);
  virtual ~ChannelingStackingMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
    
  ChannelingStackingAction*     fStackAction;
    
  G4UIcmdWithABool*   fKillCmd;
  G4UIcmdWithAString* fKCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
