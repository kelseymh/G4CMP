//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/A01/src/DetectorConstructionMessenger.cc
/// \brief Implementation of the DetectorConstructionMessenger class
//
// $Id$
// --------------------------------------------------------------
//
#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* mpga)
:fTarget(mpga)
{
    fMydetDirectory = new G4UIdirectory("/mydet/");
    fMySSDDirectory = new G4UIdirectory("/mydet/ssd/");
    fMySCIDirectory = new G4UIdirectory("/mydet/sci/");
    fMyRBSDirectory = new G4UIdirectory("/mydet/rbs/");
    fMyXtalDirectory = new G4UIdirectory("/mydet/xtal/");
    
    fMydetDirectory->SetGuidance("Channeling experiment detector setup control commands.");
    
    fSSD0XtalDistanceCmd = new G4UIcmdWithADoubleAndUnit("/mydet/ssd/setDist0X",this);
    fSSD0XtalDistanceCmd->SetGuidance("Distance Xtal and SD0.");
    fSSD0XtalDistanceCmd->SetParameterName("distSSD0Xtal",true);
    //fSSD0XtalDistanceCmd->SetRange("distSSD0Xtal>-150. && distSSD0Xtal<-25.");
    fSSD0XtalDistanceCmd->SetDefaultValue(-105.);
    fSSD0XtalDistanceCmd->SetDefaultUnit("cm");
    
    fSSD1XtalDistanceCmd = new G4UIcmdWithADoubleAndUnit("/mydet/ssd/setDist1X",this);
    fSSD1XtalDistanceCmd->SetGuidance("Distance Xtal and SD1.");
    fSSD1XtalDistanceCmd->SetParameterName("distSSD1Xtal",true);
    //fSSD1XtalDistanceCmd->SetRange("distSSD1Xtal>-20. && distSSD1Xtal<5.");
    fSSD1XtalDistanceCmd->SetDefaultValue(-5.);
    fSSD1XtalDistanceCmd->SetDefaultUnit("cm");
    
    fSSD2XtalDistanceCmd = new G4UIcmdWithADoubleAndUnit("/mydet/ssd/setDist2X",this);
    fSSD2XtalDistanceCmd->SetGuidance("Distance Xtal and SD2.");
    fSSD2XtalDistanceCmd->SetParameterName("distSSD2Xtal",true);
    //fSSD2XtalDistanceCmd->SetRange("distSSD2Xtal>-5. && distSSD2Xtal<150.");
    fSSD2XtalDistanceCmd->SetDefaultValue(+95.);
    fSSD2XtalDistanceCmd->SetDefaultUnit("cm");
    
    
    fSCIRelativeDistanceCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sci/setDist",this);
    fSCIRelativeDistanceCmd->SetGuidance("Distance between the two scintillators.");
    fSCIRelativeDistanceCmd->SetParameterName("distSCI",true);
    //fSCIRelativeDistanceCmd->SetRange("distSCI>0. && distSCI<=50.");
    fSCIRelativeDistanceCmd->SetDefaultValue(1.);
    fSCIRelativeDistanceCmd->SetDefaultUnit("cm");
    
    fSCIXtalDistanceCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sci/setDistX",this);
    fSCIXtalDistanceCmd->SetGuidance("Distance between the scintillators and the xtal");
    fSCIXtalDistanceCmd->SetParameterName("distSCIXtal",true);
    //fSCIXtalDistanceCmd->SetRange("distSCIXtal>5. && distSCIXtal<=100.");
    fSCIXtalDistanceCmd->SetDefaultValue(60.);
    fSCIXtalDistanceCmd->SetDefaultUnit("cm");
    
    
    fRBSDistanceRCmd = new G4UIcmdWithADoubleAndUnit("/mydet/rbs/setDistX",this);
    fRBSDistanceRCmd->SetGuidance("Distance between the RBS detector and the xtal.");
    fRBSDistanceRCmd->SetParameterName("distRBS",true);
    fRBSDistanceRCmd->SetRange("distRBS>1. && distRBS<=30.");
    fRBSDistanceRCmd->SetDefaultValue(20.);
    fRBSDistanceRCmd->SetDefaultUnit("cm");
    
    fRBSAngleThetaCmd = new G4UIcmdWithADoubleAndUnit("/mydet/rbs/setAngTheta",this);
    fRBSAngleThetaCmd->SetGuidance("Polar angle of the RBS detector.");
    fRBSAngleThetaCmd->SetParameterName("angRBStheta",true);
    fRBSAngleThetaCmd->SetRange("angRBStheta>=0. && angRBStheta<180.");
    fRBSAngleThetaCmd->SetDefaultValue(0.);
    fRBSAngleThetaCmd->SetDefaultUnit("deg");
    
    fRBSAnglePhiCmd = new G4UIcmdWithADoubleAndUnit("/mydet/rbs/setAngPhi",this);
    fRBSAnglePhiCmd->SetGuidance("Azimuthal angle of the RBS detector.");
    fRBSAnglePhiCmd->SetParameterName("angRBSphi",true);
    fRBSAnglePhiCmd->SetRange("angRBSphi>=0. && angRBSphi<360.");
    fRBSAnglePhiCmd->SetDefaultValue(-135.);
    fRBSAnglePhiCmd->SetDefaultUnit("deg");
    
    
    fXtalMaterialCmd = new G4UIcmdWithAString("/mydet/xtal/setMaterial",this);
    fXtalMaterialCmd->SetGuidance("Set Xtal material.");
    fXtalMaterialCmd->SetParameterName("xMat",true);
    fXtalMaterialCmd->SetDefaultValue("G4_Si");
    
    fXtalCurvatureRadiusCmd = new G4UIcmdWithADoubleAndUnit("/mydet/xtal/setCurvRadius",this);
    fXtalCurvatureRadiusCmd->SetGuidance("Set Xtal curvature radius.");
    fXtalCurvatureRadiusCmd->SetParameterName("xtalcr",true);
    fXtalCurvatureRadiusCmd->SetDefaultValue(0.);
    fXtalCurvatureRadiusCmd->SetDefaultUnit("m");

    fXtalSizeCmd = new G4UIcmdWith3VectorAndUnit("/mydet/xtal/setSize",this);
    fXtalSizeCmd->SetGuidance("Set Xtal size.");
    fXtalSizeCmd->SetParameterName("xtalSizeX","xtalSizeY","xtalSizeZ",true);
    fXtalSizeCmd->SetDefaultValue(G4ThreeVector(6.,2.,6.));
    fXtalSizeCmd->SetDefaultUnit("mm");
    
    fXtalAngleCmd = new G4UIcmdWith3VectorAndUnit("/mydet/xtal/setAngle",this);
    fXtalAngleCmd->SetGuidance("Set Xtal angles in the world.");
    fXtalAngleCmd->SetParameterName("xtalAngleX","xtalAngleY","xtalAngleZ",true);
    fXtalAngleCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalAngleCmd->SetDefaultUnit("rad");
    
    fXtalCellSizeCmd = new G4UIcmdWith3VectorAndUnit("/mydet/xtal/setCellSize",this);
    fXtalCellSizeCmd->SetGuidance("Set Xtal cell size.");
    fXtalCellSizeCmd->SetParameterName("xtalCellSizeX","xtalCellSizeY","xtalCellSizeZ",true);
    fXtalCellSizeCmd->SetDefaultValue(G4ThreeVector(5.431,5.431,5.431));
    fXtalCellSizeCmd->SetDefaultUnit("angstrom");
    
    fXtalCellAngleCmd = new G4UIcmdWith3VectorAndUnit("/mydet/xtal/setCellAngle",this);
    fXtalCellAngleCmd->SetGuidance("Set Xtal cell angles in the world.");
    fXtalCellAngleCmd->SetParameterName("xtalCellAngleX","xtalCellAngleY","xtalCellAngleZ",true);
    fXtalCellAngleCmd->SetDefaultValue(G4ThreeVector(90.,90.,90.));
    fXtalCellAngleCmd->SetDefaultUnit("deg");

    fAddRBSDetector = new G4UIcmdWithAString("/mydet/addrbsd",this);
    fAddRBSDetector->SetGuidance("Add RBS detector.");
    fAddRBSDetector->SetParameterName("addrbsd",true);

    fAddScintillator = new G4UIcmdWithAString("/mydet/addsci",this);
    fAddScintillator->SetGuidance("Add scintillators.");
    fAddScintillator->SetParameterName("addsci",true);

    fAddSiliconDetector = new G4UIcmdWithAString("/mydet/addsd",this);
    fAddSiliconDetector->SetGuidance("Add SD detectors.");
    fAddSiliconDetector->SetParameterName("addsd",true);

    fAddXtalTarget = new G4UIcmdWithAString("/mydet/addxtal",this);
    fAddXtalTarget->SetGuidance("Add Xtal.");
    fAddXtalTarget->SetParameterName("addxtal",true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
    delete fMydetDirectory;
    
    delete fSSD0XtalDistanceCmd;
    delete fSSD1XtalDistanceCmd;
    delete fSSD2XtalDistanceCmd;
    
    delete fSCIRelativeDistanceCmd;
    delete fSCIXtalDistanceCmd;
    
    delete fRBSDistanceRCmd;
    delete fRBSAngleThetaCmd;
    delete fRBSAnglePhiCmd;
    
    delete fXtalMaterialCmd;
    delete fXtalCurvatureRadiusCmd;
    delete fXtalSizeCmd;
    delete fXtalAngleCmd;
    delete fXtalCellSizeCmd;
    delete fXtalCellAngleCmd;
    
    delete fAddRBSDetector;
    delete fAddScintillator;
    delete fAddSiliconDetector;
    delete fAddXtalTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
    if(command==fAddRBSDetector ){
        fTarget->AddRBSDetector();
    }
    if(command==fAddScintillator ){
        fTarget->AddScintillators();
    }
    if(command==fAddSiliconDetector ){
        fTarget->AddSiliconStripDetectors();
    }
    if(command==fAddXtalTarget ){
        fTarget->AddXtalTarget();
    }

    
    if(command==fSSD0XtalDistanceCmd ){
        fTarget->SetSSD0XtalDistance(fSSD0XtalDistanceCmd->GetNewDoubleValue(newValue));
    }
    if(command==fSSD0XtalDistanceCmd ){
        fTarget->SetSSD1XtalDistance(fSSD1XtalDistanceCmd->GetNewDoubleValue(newValue));
    }
    if(command==fSSD0XtalDistanceCmd ){
        fTarget->SetSSD2XtalDistance(fSSD1XtalDistanceCmd->GetNewDoubleValue(newValue));
    }
    
    
    if(command==fSCIRelativeDistanceCmd ){
        fTarget->SetSCIRelativeDistance(fSCIRelativeDistanceCmd->GetNewDoubleValue(newValue));
    }
    if(command==fSCIXtalDistanceCmd ){
        fTarget->SetSCIXtalDistance(fSCIXtalDistanceCmd->GetNewDoubleValue(newValue));
    }
    
    
    if(command==fRBSDistanceRCmd ){
        fTarget->SetRBSDistanceR(fRBSDistanceRCmd->GetNewDoubleValue(newValue));
    }
    if(command==fRBSAngleThetaCmd ){
        fTarget->SetRBSAngleTheta(fRBSAngleThetaCmd->GetNewDoubleValue(newValue));
    }
    if(command==fRBSAnglePhiCmd ){
        fTarget->SetRBSAnglePhi(fRBSAnglePhiCmd->GetNewDoubleValue(newValue));
    }
    
    
    if(command==fXtalMaterialCmd ){
        fTarget->SetXtalMaterial(newValue);
    }
    if(command==fXtalCurvatureRadiusCmd ){
        fTarget->SetXtalCurvatureRadius(fXtalCurvatureRadiusCmd->GetNewDoubleValue(newValue));
    }
    if(command==fXtalSizeCmd ){
        fTarget->SetXtalSize(fXtalSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalAngleCmd ){
        fTarget->SetXtalAngle(fXtalAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellSizeCmd ){
        fTarget->SetXtalCellSize(fXtalCellSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellAngleCmd ){
        fTarget->SetXtalCellAngle(fXtalCellAngleCmd->GetNew3VectorValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String DetectorConstructionMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    
    if( command==fSSD0XtalDistanceCmd ){
        cv = fSSD0XtalDistanceCmd->ConvertToString(fTarget->GetSSD0XtalDistance(),"cm");
    }
    if( command==fSSD1XtalDistanceCmd ){
        cv = fSSD1XtalDistanceCmd->ConvertToString(fTarget->GetSSD1XtalDistance(),"cm");
    }
    if( command==fSSD1XtalDistanceCmd ){
        cv = fSSD1XtalDistanceCmd->ConvertToString(fTarget->GetSSD1XtalDistance(),"cm");
    }
    
    
    if( command==fSCIRelativeDistanceCmd ){
        cv = fSCIRelativeDistanceCmd->ConvertToString(fTarget->GetSCIRelativeDistance(),"cm");
    }
    if( command==fSCIXtalDistanceCmd ){
        cv = fSCIXtalDistanceCmd->ConvertToString(fTarget->GetSCIXtalDistance(),"cm");
    }
    
    
    if( command==fRBSDistanceRCmd ){
        cv = fRBSDistanceRCmd->ConvertToString(fTarget->GetRBSDistanceR(),"cm");
    }
    if( command==fRBSAngleThetaCmd ){
        cv = fRBSAngleThetaCmd->ConvertToString(fTarget->GetRBSAngleTheta(),"deg");
    }
    if( command==fRBSAnglePhiCmd ){
        cv = fRBSAnglePhiCmd->ConvertToString(fTarget->GetRBSAnglePhi(),"deg");
    }
    
    
    if( command==fXtalMaterialCmd ){
        cv = fTarget->GetXtalMaterial();
    }
    if( command==fXtalCurvatureRadiusCmd ){
        cv = fXtalCurvatureRadiusCmd->ConvertToString(fTarget->GetXtalCurvatureRadius(),"m");
    }
    if( command==fXtalSizeCmd ){
        cv = fXtalSizeCmd->ConvertToString(fTarget->GetXtalSize(),"mm");
    }
    if( command==fXtalAngleCmd ){
        cv = fXtalAngleCmd->ConvertToString(fTarget->GetXtalAngle(),"rad");
    }
    if( command==fXtalCellSizeCmd ){
        cv = fXtalCellSizeCmd->ConvertToString(fTarget->GetXtalCellSize(),"angstrom");
    }
    if( command==fXtalCellAngleCmd ){
        cv = fXtalCellAngleCmd->ConvertToString(fTarget->GetXtalCellAngle(),"deg");
    }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
