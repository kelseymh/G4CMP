/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingDetectorMessenger.hh
/// \brief Definition of the ChannelingDetectorMessenger class
//
// $Id: 96d0a882244d54eb004c6f0b35650603fee30991 $
// --------------------------------------------------------------
//
#ifndef ChannelingDetectorMessenger_h
#define ChannelingDetectorMessenger_h 1

class ChannelingDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"

class ChannelingDetectorMessenger: public G4UImessenger
{
  public:
    ChannelingDetectorMessenger(ChannelingDetectorConstruction* mpga);
    ~ChannelingDetectorMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    ChannelingDetectorConstruction * fTarget;

    G4UIdirectory* fMydetDirectory;
    G4UIdirectory* fMySSDDirectory;
    G4UIdirectory* fMySCIDirectory;
    G4UIdirectory* fMyRBSDirectory;
    G4UIdirectory* fMyXtalDirectory;

    G4UIcmdWithAString*  fWorldMaterialCmd;

    G4UIcmdWithAString*  fAddRBSDetector;
    G4UIcmdWithAString*  fAddScintillator;
    G4UIcmdWithAString*  fAddSiliconDetector;
    G4UIcmdWithAString*  fAddXtalTarget;

    G4UIcmdWith3VectorAndUnit*  fSSDSizeCmd;
    G4UIcmdWithADoubleAndUnit*  fSSD0XtalDistanceCmd;
    G4UIcmdWithADoubleAndUnit*  fSSD1XtalDistanceCmd;
    G4UIcmdWithADoubleAndUnit*  fSSD2XtalDistanceCmd;
    
    G4UIcmdWith3VectorAndUnit*  fSCISizeCmd;
    G4UIcmdWithADoubleAndUnit*  fSCIRelativeDistanceCmd;
    G4UIcmdWithADoubleAndUnit*  fSCIXtalDistanceCmd;
    
    G4UIcmdWithADoubleAndUnit*  fRBSDistanceRCmd;
    G4UIcmdWithADoubleAndUnit*  fRBSAngleThetaCmd;
    G4UIcmdWithADoubleAndUnit*  fRBSAnglePhiCmd;
    
    G4UIcmdWithAString*  fXtalMaterialCmd;
    G4UIcmdWith3VectorAndUnit* fXtalCurvatureRadiusCmd;
    G4UIcmdWith3VectorAndUnit* fXtalSizeCmd;
    G4UIcmdWith3VectorAndUnit* fXtalAngleCmd;
    G4UIcmdWith3VectorAndUnit* fXtalCellSizeCmd;
    G4UIcmdWith3VectorAndUnit* fXtalCellAngleCmd;
    G4UIcmdWithADoubleAndUnit* fXtalCellThermalVibration;
};

#endif


