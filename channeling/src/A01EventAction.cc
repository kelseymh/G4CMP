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
    fDHC_ID = -1;

    fVerboseLevel = 0;
    fMessenger = new A01EventActionMessenger(this);
}

A01EventAction::~A01EventAction()
{
    delete fMessenger;
}

void A01EventAction::BeginOfEventAction(const G4Event* evt){
}

void A01EventAction::EndOfEventAction(const G4Event* evt)
{
    
//    if(fDHC_ID==-1) {
//        G4String colName;
//        G4SDManager* SDman = G4SDManager::GetSDMpointer();
//        fDHC_ID = SDman->GetCollectionID(colName="telescope/telescopeColl");
//    }
//    
//    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
//    A01DriftChamberHitsCollection* fDHC = 0;
//    
//    if(HCE)
//    {
//        G4VHitsCollection* aHC = HCE->GetHC(fDHC_ID);
//        //fDHC = dynamic_cast<A01DriftChamberHitsCollection*>(aHC);
//        fDHC = (A01DriftChamberHitsCollection*)(aHC);
//    }
//    
//    // Diagnostics
//    
//    if (fVerboseLevel==0 || evt->GetEventID() % fVerboseLevel != 0) return;
//    
//    G4PrimaryParticle* primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
//    G4cout << G4endl
//    << ">>> Event " << evt->GetEventID() << " >>> Simulation truth : "
//    << primary->GetG4code()->GetParticleName()
//    << " " << primary->GetMomentum() << G4endl;
//    
//    
//    {
//        if(fDHC)
//        {
//            int n_hit = fDHC->entries();
//            G4cout << "Drift Chamber has " << n_hit << " hits." << G4endl;
//            for(int i2=0;i2<5;i2++)
//            {
//                for(int i1=0;i1<n_hit;i1++)
//                {
//                    A01DriftChamberHit* aHit = (*fDHC)[i1];
//                    if(aHit->GetLayerID()==i2) aHit->Print();
//                }
//            }
//        }
//    }
}


