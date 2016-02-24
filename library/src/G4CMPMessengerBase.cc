/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
////////////////////////////////////////////////////////////////////////
//  File:        G4CMPMessengerBase.cc                                 //     
//  Description: Base class for all G4CMP GEANT4 Messengers. Provides //
//		 common functionality and some diagnostic utilities.  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 August 2014                                        //
//////////////////////////////////////////////////////////////////////// 

#include "G4CMPMessengerBase.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"


// Constructor looks to see if directory exists before creating it

G4CMPMessengerBase::G4CMPMessengerBase(const G4String& path,
				       const G4String& desc)
  : localCmdDir(false), cmdDir(0) {
  CreateDirectory(path, desc);
}

// Destructor only deletes directory if it was created here

G4CMPMessengerBase::~G4CMPMessengerBase() {
  if (localCmdDir) delete cmdDir;
}


// Create new command directory, or reference existing version

void G4CMPMessengerBase::CreateDirectory(const char* path, const char* desc) {
  // If no UI created, everything gets ignored
  G4UImanager* UIman = G4UImanager::GetUIpointer();
  if (!UIman) return;

  // Prepend /g4cmp/ if user specified "relative path" (e.g., module name)
  G4String fullPath = path;
  if (fullPath(0) != '/') fullPath.prepend("/g4cmp/");
  if (fullPath(fullPath.length()-1) != '/') fullPath.append("/");

  // See if input path has already been registered
  G4UIcommand* foundPath = UIman->GetTree()->FindPath(fullPath);
  if (foundPath) cmdDir = dynamic_cast<G4UIdirectory*>(foundPath);

  if (!cmdDir) {		// Create local deletable directory
    localCmdDir = true;
    cmdDir = new G4UIdirectory(fullPath);
    cmdDir->SetGuidance(desc);
  }
}
