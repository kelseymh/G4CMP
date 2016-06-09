/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef G4CMPChannelingSlow_h
#define G4CMPChannelingSlow_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "XLatticeManager3.hh"
#include "XVCrystalCharacteristic.hh"
#include "XCrystalIntegratedDensityHub.hh"

#include "G4VPhysicalVolume.hh"
#include "XPhysicalLattice.hh"
#include "ChannelingParticleUserInfo.hh"

class G4CMPChannelingSlow : public G4VDiscreteProcess
{
public:
    
    G4CMPChannelingSlow(const G4String& processName = "G4CMPChannelingSlow" );
    
    virtual ~G4CMPChannelingSlow();
    
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
        
    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    
    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    
protected:
    
    virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* );
    
private:
    
    G4double GetChannelingMeanFreePath(const G4Track&);
    
    G4bool IsUnderCoherentEffect(const G4Track&);

public:
    XVCrystalCharacteristic* GetPotential();
    void SetPotential(XVCrystalCharacteristic*);

    XVCrystalCharacteristic* GetElectricField();
    void SetElectricField(XVCrystalCharacteristic*);

    XCrystalIntegratedDensityHub* GetIntegratedDensity();
    void SetIntegratedDensity(XCrystalIntegratedDensityHub*);

    void SetFileName(const G4String&);
    G4String GetFileName();
    void PrintChannelingPropertiesOnFile(const G4Track&);
    
private:
    void UpdateParameters(const G4Track&);
    void UpdatePosition(const G4Track&);
    void UpdateMomentum(const G4Track&);
    void ResetDensity(const G4Track&);

    G4double ComputeCriticalEnergyBent(const G4Track&);
    G4double ComputeCriticalEnergyMinimumBent(const G4Track&);
    G4double ComputePotentialEnergyBent(const G4Track&);
    G4ThreeVector ComputeTransverseEnergyBent(const G4Track&);
    
    G4ThreeVector ComputePositionInTheCrystal(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForVolume(G4StepPoint*,const G4Track&);
    G4StepPoint* CheckStepPointLatticeForPosition(G4StepPoint*,const G4Track&);
    
    G4ThreeVector ComputeTransverseEnergy(const G4Track&);
    G4ThreeVector ComputeKineticEnergy(const G4Track&);
    G4ThreeVector ComputePotentialEnergy(const G4Track&);
    G4ThreeVector ComputeCentrifugalEnergy(const G4Track&,G4ThreeVector);
    G4ThreeVector ComputeMomentum(const G4Track&,G4StepPoint*);
    
    G4double ComputeCriticalEnergyMaximum(const G4Track&);
    G4double ComputeCriticalEnergyMinimum(const G4Track&);
    G4double ComputeCriticalAngle(const G4Track&);
    G4double ComputeCriticalRadius(const G4Track&);
    G4double ComputeOscillationPeriod(const G4Track&);
    G4ThreeVector ComputePotentialWellCentre(const G4Track&);
    
    G4bool HasLattice(const G4Track&);
    G4bool HasLatticeOnBoundary(const G4Track&);
    G4bool ParticleIsNegative(const G4Track&);
    G4bool ParticleIsCrossingPlane(const G4Track&);
    G4bool ParticleIsNotOnBoundary(const G4Track&);

    void ComputeCrystalCharacteristic(const G4Track&);
    void PrintCrystalCharacteristicsOnFiles(const G4Track&);

private:
    //binding methods
    XPhysicalLattice* GetXPhysicalLattice(const G4Track&);
    G4VPhysicalVolume* GetVolume(const G4Track&);
    ChannelingParticleUserInfo* GetInfo(const G4Track&);
    G4ParticleDefinition* GetParticleDefinition(const G4Track& aTrack);
    
private:
    // hide assignment operator as private
    G4CMPChannelingSlow(G4CMPChannelingSlow&);
    G4CMPChannelingSlow& operator=(const G4CMPChannelingSlow& right);
    
private:
    XLatticeManager3* fLatticeManager;
      
    XVCrystalCharacteristic* fPotentialEnergy;
    XVCrystalCharacteristic* fElectricField;
    XCrystalIntegratedDensityHub* fIntegratedDensity;

    std::ofstream fFileOut;
    G4String fFileName;
};

#endif










