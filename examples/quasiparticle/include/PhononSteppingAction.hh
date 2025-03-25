/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef PhononSteppingAction_hh
#define PhononSteppingAction_hh 1

#include "G4UserSteppingAction.hh"

#include <fstream>

class G4Step;

class PhononSteppingAction : public G4UserSteppingAction
{
public:

  PhononSteppingAction();
  virtual ~PhononSteppingAction();
  virtual void UserSteppingAction(const G4Step* step);
  void ExportStepInformation( const G4Step * step );
  
private:

  //Step info output file
  std::ofstream fOutputFile;
  
  
  
  
};

#endif

