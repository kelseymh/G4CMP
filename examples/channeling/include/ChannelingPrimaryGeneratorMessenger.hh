/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingPrimaryGeneratorMessenger.hh
/// \brief Definition of the ChannelingPrimaryGeneratorMessenger class
//
// $Id: 18dca17fa9965f8c8555c3fd48e30ced11f6c0ae $
// --------------------------------------------------------------
//
#ifndef ChannelingPrimaryGeneratorMessenger_h
#define ChannelingPrimaryGeneratorMessenger_h 1

class ChannelingPrimaryGeneratorAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIdirectory;

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

class ChannelingPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    ChannelingPrimaryGeneratorMessenger(ChannelingPrimaryGeneratorAction* mpga);
    virtual ~ChannelingPrimaryGeneratorMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    ChannelingPrimaryGeneratorAction* fTarget;

    G4UIcmdWithAString* fDivergenceDistribution;
    G4UIcmdWithADoubleAndUnit* fCutX;
    G4UIcmdWithADoubleAndUnit* fCutY;
    G4UIcmdWithADoubleAndUnit* fDivergenceX;
    G4UIcmdWithADoubleAndUnit* fDivergenceY;
    G4UIcmdWithADoubleAndUnit* fWidthX;
    G4UIcmdWithADoubleAndUnit* fWidthY;

};

#endif


