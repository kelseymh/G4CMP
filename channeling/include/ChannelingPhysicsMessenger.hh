/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingPhysicsMessenger.hh
/// \brief Definition of the ChannelingPhysicsMessenger class
//
// $Id: fcf6ba4fe2350a27cdf5750731603db736ed51b4 $
// --------------------------------------------------------------
//
#ifndef ChannelingPhysicsMessenger_h
#define ChannelingPhysicsMessenger_h 1

class ChannelingPhysicsList;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

#include "G4UImessenger.hh"
#include "globals.hh"

class ChannelingPhysicsMessenger: public G4UImessenger
{
  public:
    ChannelingPhysicsMessenger(ChannelingPhysicsList* mpga);
    virtual ~ChannelingPhysicsMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    ChannelingPhysicsList* fTarget;

    G4UIcmdWithAString* fFileNameCmd;
    G4UIcmdWithAString* fScatteringType;

    G4UIcmdWithABool* fChannelingCmd;
    G4UIcmdWithABool* fWrapperCmd;
    G4UIcmdWithABool* fDecayCmd;

    G4UIdirectory* fMyDirectory;
};

#endif


