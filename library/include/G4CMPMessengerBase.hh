/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPMessengerBase_hh
#define G4CMPMessengerBase_hh 1
// $Id$
////////////////////////////////////////////////////////////////////////
//  File:        G4CMPMessengerBase.hh                                //     
//  Description: Base class for all G4CMP GEANT4 Messengers. Provides //
//		 common functionality and some diagnostic utilities.  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 August 2014                                        //
//////////////////////////////////////////////////////////////////////// 

#include "G4UImessenger.hh"
#include "globals.hh"
#include <iosfwd>

class G4UIdirectory;
class G4UIcmdWithAnInteger;


class G4CMPMessengerBase : public G4UImessenger {
public:
  G4CMPMessengerBase(const G4String& path, const G4String& desc);
  virtual ~G4CMPMessengerBase();

  // NOTE:  Subclasses must define this command needed by G4UImanager
  virtual void SetNewValue(G4UIcommand* command, G4String newValue) = 0;

protected:
  // Create G4UIcommand (or subclass) within current command path
  template <class T>
  T* CreateCommand(const G4String& cmd, const G4String& desc) const;

private:
  // Configuration function used by constructor
  void CreateDirectory(const char* path, const char* desc);

private:
  G4bool localCmdDir;		// Flag if directory was created or found
  G4UIdirectory* cmdDir;
};

#include "G4CMPMessengerBase.icc"

#endif /* G4CMPMessengerBase_hh */
