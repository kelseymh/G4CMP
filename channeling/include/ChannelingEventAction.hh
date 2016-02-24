/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingEventAction.hh
/// \brief Definition of the ChannelingEventAction class
//
// $Id: 938ad2a7905e199e200a3acd9e372bfcb9bb1eda $
// --------------------------------------------------------------
//
#ifndef ChannelingEventAction_h
#define ChannelingEventAction_h 1


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

class ChannelingEventMessenger;

class ChannelingEventAction : public G4UserEventAction
{
  public:
    ChannelingEventAction();
    virtual ~ChannelingEventAction();

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

    ChannelingEventMessenger* fMessenger;
    G4int fVerboseLevel;

  public:
    inline void SetVerbose(G4int val) { fVerboseLevel = val; }
    inline G4int GetVerbose() const { return fVerboseLevel; }
};

#endif
