/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPPhysicsList.cc
/// \brief Definition of the G4CMPPhysicsList class
//
// $Id$
//
// 20140328  Change to ModularPhysicsList, use new G4CMPPhysics for processes
// 20170525  Destructor should be virtual, add default copy/move semantics

#ifndef G4CMPPhysicsList_h
#define G4CMPPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4CMPPhysicsList: public G4VModularPhysicsList {
public:
  G4CMPPhysicsList(G4int verbose=0);
  virtual ~G4CMPPhysicsList();

  // Use default copy/move semantics
  G4CMPPhysicsList(const G4CMPPhysicsList&) = default;
  G4CMPPhysicsList(G4CMPPhysicsList&&) = default;
  G4CMPPhysicsList& operator=(const G4CMPPhysicsList&) = default;
  G4CMPPhysicsList& operator=(G4CMPPhysicsList&&) = default;

public:
  void SetCuts();
};

#endif
