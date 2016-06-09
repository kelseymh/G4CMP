/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file electromagnetic/TestEm5/include/ChannelingStepMessenger.hh
/// \brief Definition of the ChannelingStepMessenger class
//
// $Id: 19df442e5abb311b15541dae39e7993f0f2d23d4 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ChannelingStepMessenger_h
#define ChannelingStepMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class ChannelingStepLimiter;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChannelingStepMessenger: public G4UImessenger
{
  public:
    ChannelingStepMessenger(ChannelingStepLimiter*);
   ~ChannelingStepMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ChannelingStepLimiter* fChannelingStepLimiter;
    G4UIcmdWithADoubleAndUnit* fChannelingStepLimiterCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
