/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingEventMessenger.hh
/// \brief Definition of the ChannelingEventMessenger class
//
// $Id: faa7feacf4e91f8c5cb94284de8654fcf98a8153 $
// --------------------------------------------------------------
//
#ifndef ChannelingEventMessenger_h
#define ChannelingEventMessenger_h 1

class ChannelingEventAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIdirectory;

#include "G4UImessenger.hh"
#include "globals.hh"

class ChannelingEventMessenger: public G4UImessenger
{
  public:
    ChannelingEventMessenger(ChannelingEventAction* mpga);
    virtual ~ChannelingEventMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    ChannelingEventAction* fTarget;

    G4UIdirectory* fMyDirectory;
    G4UIcmdWithAnInteger* fVerboseCmd;
    G4UIcmdWithAString* fFileNameCmd;

};

#endif


