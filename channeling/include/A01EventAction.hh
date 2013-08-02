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
/// \file analysis/A01/include/A01EventAction.hh
/// \brief Definition of the A01EventAction class
//
// $Id$
// --------------------------------------------------------------
//
#ifndef A01EventAction_h
#define A01EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>

#ifdef ROOT
#include "TFile.h"
#include "TTree.h"
struct ROOT_save
{
    Double_t x;
    Double_t y;
    Double_t z;
    Double_t en;
    UInt_t id;
    UInt_t layer;
};
#endif

class A01EventActionMessenger;

class A01EventAction : public G4UserEventAction
{
  public:
    A01EventAction();
    virtual ~A01EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void SetFileName(const G4String&);
    G4String GetFileName();
    
  private:
    G4int fSD_ID;
    G4int fSCI_ID;
    G4int fRBSD_ID;

    G4String fFileName;
    std::ofstream fFileOutSCI;
    std::ofstream fFileOutSD;
    std::ofstream fFileOutRBSD;

#ifdef ROOT
    TFile *fRootFile;
    TTree *fTree;
    const ROOT_save fDetectorSave;
#endif

    A01EventActionMessenger* fMessenger;
    G4int fVerboseLevel;

  public:
    inline void SetVerbose(G4int val) { fVerboseLevel = val; }
    inline G4int GetVerbose() const { return fVerboseLevel; }
};

#endif
