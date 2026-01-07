 /***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file ValidationSteppingAction.hh
/// \brief Definition of the class governing what we do at each step in the
///  validation example

#ifndef ValidationSteppingAction_hh
#define ValidationSteppingAction_hh 1

#include "G4UserSteppingAction.hh"

#include <fstream>

class G4Step;

class ValidationSteppingAction : public G4UserSteppingAction
{
public:

  ValidationSteppingAction();
  virtual ~ValidationSteppingAction();
  virtual void UserSteppingAction(const G4Step* step);
  void ExportStepInformation( const G4Step * step );
  
private:

  //Step info output file
  std::ofstream fOutputFile;
  
  
  
  
};

#endif

