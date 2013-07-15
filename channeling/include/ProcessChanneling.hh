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
//

#ifndef ProcessChanneling_h
#define ProcessChanneling_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "XLatticeManager3.hh"
#include "XVCrystalCharacteristic.hh"
#include "XCrystalIntegratedDensityHub.hh"

#include "G4VPhysicalVolume.hh"
#include "XPhysicalLattice.hh"
#include "ChannelingParticleUserInfo.hh"

class ProcessChanneling : public G4VDiscreteProcess
{
public:
    
    ProcessChanneling(const G4String& processName = "ProcessChanneling" );
    
    virtual ~ProcessChanneling();
    
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
    
    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    
    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    
protected:
    virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* );
    
public:
    XVCrystalCharacteristic* GetPotential();
    void SetPotential(XVCrystalCharacteristic*);
    
    XCrystalIntegratedDensityHub* GetIntegratedDensity();
    void SetIntegratedDensity(XCrystalIntegratedDensityHub*);

    void SetFileName(const G4String&);
    G4String GetFileName();
    void PrintChannelingPropertiesOnFile(const G4Track&);
    
private:
    G4double GetChannelingMeanFreePath(const G4Track&);
    
    G4bool IsUnderCoherentEffect(const G4Track&);

    G4bool HasLattice(const G4Track&);
    G4bool HasLatticeOnBoundary(const G4Track&);

    void UpdatePositionMomentumDensity(const G4Track&);
    void UpdatePosition(const G4Track&);
    void UpdateMomentum(const G4Track&);
    void UpdateDensity(const G4Track&);

    G4ThreeVector ComputeTransverseEnergy(const G4Track&);
    G4ThreeVector ComputeTransverseEnergyBent(const G4Track&);
    
    G4ThreeVector ComputeChannelingOutgoingMomentum(const G4Track&);
    G4ThreeVector ComputeVolumeReflectionOutgoingMomentum(const G4Track&);

    G4ThreeVector ComputePostStepPositionInTheChannel(const G4Track&);
    G4ThreeVector ComputePositionInTheCrystal(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForVolume(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForPosition(G4StepPoint*,const G4Track&);

    void ComputeCrystalCharacteristic(const G4Track&);
    void PrintCrystalCharacteristicsOnFiles(const G4Track&);
    
    G4double ComputeCriticalEnergy(const G4Track&);
    G4double ComputeCriticalEnergyBent(const G4Track&);
    G4double ComputeCriticalAngle(const G4Track&);
    G4double ComputeOscillationPeriod(const G4Track&);
    
private:
    //binding methods
    XPhysicalLattice* GetXPhysicalLattice(const G4Track&);
    G4VPhysicalVolume* GetVolume(const G4Track&);
    ChannelingParticleUserInfo* GetInfo(const G4Track&);
    G4ParticleDefinition* GetParticleDefinition(const G4Track& aTrack);
    
private:
    // hide assignment operator as private
    ProcessChanneling(ProcessChanneling&);
    ProcessChanneling& operator=(const ProcessChanneling& right);
    
private:
    XLatticeManager3* fLatticeManager;
      
    XVCrystalCharacteristic* fPotentialEnergy;
    XCrystalIntegratedDensityHub* fIntegratedDensity;

    std::ofstream fFileOut;
    G4String fFileName;
};

#endif










