/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef CHARGEFIELDMESSENGER_HH
#define CHARGEFIELDMESSENGER_HH 1

class G4UIdirectory;
class G4UIcommandTree;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWithoutParameter;
class ChargeEMField;

#include "G4UImessenger.hh"
#include "globals.hh"

class ChargeFieldMessenger: public G4UImessenger
{
  public:
    ChargeFieldMessenger(ChargeEMField* field);
    ~ChargeFieldMessenger();
    void SetNewValue(G4UIcommand* command, G4String NewValue);
  private:
    ChargeEMField* ChargeField;
    G4UIdirectory* fieldDir;
    G4UIcommandTree* fieldTree;
    G4UIcmdWithoutParameter* buildcmd;
    G4UIcmdWithABool* constFieldcmd;
    G4UIcmdWithADoubleAndUnit* constFieldMagcmd;
    G4UIcmdWith3Vector* constFieldDircmd;
    G4UIcmdWithAString* fieldFileNamecmd;
    G4bool localFieldDir;
  /*protected:
  // Create G4UIcommand (arbitrary subclass) within current command path
  template <class T>
  T* CreateCommand(const G4String& commandName, const G4String& description) {
    G4String path = fieldDir ? fieldDir->GetCommandPath() : "";
    path += commandName;

    T* theCmd = new T(path.c_str(), this);	// <T> must be G4UIcommand!
    theCmd->SetGuidance(description.c_str());
    theCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    return theCmd;
  }
  */
};

#endif
