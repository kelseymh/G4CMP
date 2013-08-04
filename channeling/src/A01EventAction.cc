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
/// \file analysis/A01/src/A01EventAction.cc
/// \brief Implementation of the A01EventAction class
//
// $Id$
// --------------------------------------------------------------
//

#include "A01EventAction.hh"
#include "A01EventActionMessenger.hh"

#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "A01DriftChamberHit.hh"

#include "TrackingAction.hh"

A01EventAction::A01EventAction()
{
    fSD_ID = -1;
    fSCI_ID = -1;
    fRBSD_ID = -1;
    
    fFileName = "";
    
    fVerboseLevel = 0;
    fMessenger = new A01EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

A01EventAction::~A01EventAction()
{
    fFileOutSCI.close();
    fFileOutSD.close();
    fFileOutRBSD.close();
    
#ifdef ROOT
    fTree->Write();
    fTree->Close();
#endif
    
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void A01EventAction::SetFileName(const G4String& vFileName){
    fFileOutSCI.close();
    fFileOutSD.close();
    fFileOutRBSD.close();
    
    fFileName = vFileName;
    
    G4String filename;
    
    filename = fFileName + "_SCI.txt";
    fFileOutSCI.open(filename.c_str());
    
    filename = fFileName + "_SD.txt";
    fFileOutSD.open(filename.c_str());
    
    filename = fFileName + "_RBSD.txt";
    fFileOutRBSD.open(filename.c_str());
    
#ifdef ROOT
    std::string vRootFileName = fFileName + ".root";
    fRootFile = new TFile(vRootFileName.c_str(),"RECREATE");
    fTree = new TTree("t","G4CMP Channeling");
    fTree->Branch("det",&fDetectorSave,"x/D:y/D:z/D:en/D:id/i:layer/i");
#endif
    
    fFileOutSCI << "hit,det,posx,posy,posz,en,id" << std::endl;
    //    fFileOutSD << "hit,det,posx,posy,posz,en,id" << std::endl;
    fFileOutSD << "x1,y1,z1,x2,y2,z2,x3,y3,z3,s1,s2,id" << std::endl;
    fFileOutRBSD << "hit,det,posx,posy,posz,en,id" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String A01EventAction::GetFileName(){
    return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void A01EventAction::BeginOfEventAction(const G4Event*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void A01EventAction::EndOfEventAction(const G4Event* evt)
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    
    if(fSD_ID==-1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            fSD_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    
    if(fSCI_ID==-1) {
        G4String sciName;
        if(SDman->FindSensitiveDetector(sciName="scintillator",0)){
            fSCI_ID = SDman->GetCollectionID(sciName="scintillator/collection");
        }
    }
    
    if(fRBSD_ID==-1) {
        G4String rbsdName;
        if(SDman->FindSensitiveDetector(rbsdName="detectorRBS",0)){
            fRBSD_ID = SDman->GetCollectionID(rbsdName="detectorRBS/collection");
        }
    }
    
    A01DriftChamberHitsCollection* fSD = 0;
    A01DriftChamberHitsCollection* fSCI = 0;
    A01DriftChamberHitsCollection* fRBSD = 0;
    
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    
    if(HCE)
    {
        //fSD = dynamic_cast<A01DriftChamberHitsCollection*>(aHC);

        if(fSCI_ID != -1){
            G4VHitsCollection* aHCSCI = HCE->GetHC(fSCI_ID);
            fSCI = (A01DriftChamberHitsCollection*)(aHCSCI);
        }
        
        if(fSD_ID != -1){
            G4VHitsCollection* aHCSD = HCE->GetHC(fSD_ID);
            fSD = (A01DriftChamberHitsCollection*)(aHCSD);
        }
        
        if(fRBSD_ID != -1){
            G4VHitsCollection* aHCRBSD = HCE->GetHC(fRBSD_ID);
            fRBSD = (A01DriftChamberHitsCollection*)(aHCRBSD);
        }
    }
    
    // Diagnostics
        
    G4PrimaryParticle* primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
    if(fmod(evt->GetEventID(),1)==0){
        G4cout << G4endl
        << ">>> Event " << evt->GetEventID() << " >>> Simulation truth : "
        << primary->GetG4code()->GetParticleName()
        << " " << primary->GetMomentumDirection() << G4endl;
    }
        
    if(fSD && fFileOutSD)
    {
        int n_hit_sd = fSD->entries();
        for(int i2=0;i2<3;i2++){
            bool bVoid = true;
            for(int i1=0;i1<n_hit_sd;i1++)
            {
                A01DriftChamberHit* aHit = (*fSD)[i1];
                if(aHit->GetLayerID()==i2 && bVoid == true) {
                    fFileOutSD << aHit->GetWorldPos().x() << "," << aHit->GetWorldPos().y() << "," << aHit->GetWorldPos().z() << ",";
                    bVoid = false;
                }
            }
            if(bVoid == true){
                fFileOutSD << "0.,0.,0.,";
            }
        }
    }

    int hitA = 0;
    int hitB = 0;

    if(fSCI && fFileOutSD)
    {
        int n_hit_sd = fSCI->entries();
            for(int i1=0;i1<n_hit_sd;i1++)
            {
                A01DriftChamberHit* aHit = (*fSCI)[i1];
                if(aHit->GetLayerID()==0) {
                    hitA++;
                }
                if(aHit->GetLayerID()==1) {
                    hitB++;
                }
            }
    }

    fFileOutSD << hitA << "," << hitB << ",";

    fFileOutSD <<  evt->GetEventID() <<  std::endl;
    
    
    
    if(fRBSD && fFileOutRBSD)
    {
        int n_hit_sd = fRBSD->entries();
        for(int i1=0;i1<n_hit_sd;i1++)
        {
            A01DriftChamberHit* aHit = (*fRBSD)[i1];
            fFileOutRBSD << i1 << "," << aHit->GetLayerID() << "," << aHit->GetWorldPos().x() << "," << aHit->GetWorldPos().y() << "," << aHit->GetWorldPos().z()<< "," << aHit->GetEnergy() << ","<< evt->GetEventID() <<  std::endl;
        }
    }
    
    
#ifdef ROOT
    if((fSCI || fSD || fRBSD)){
        fDetectorSave.x = aHit->GetWorldPos().x();
        fDetectorSave.y = aHit->GetWorldPos().y();
        fDetectorSave.z = aHit->GetWorldPos().z();
        fDetectorSave.en = aHit->GetEnergy();
        fDetectorSave.id = evt->GetEventID();
        fDetectorSave.layer = aHit->GetLayerID();
        fTree->Fill();
    }
#endif
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
