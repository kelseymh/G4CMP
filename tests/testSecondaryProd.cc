/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Usage: testSecondaryProd <trackFile> [combiningLength unit]
//
// Exercises G4CMPSecondaryProduction and G4CMPStepAccumulator on a preset
// track with hits.
//
// <trackFile>:	Name of text file describing all the steps of a track; the
// 	first line describes the track kinematics, while subesequent lines
// 	describe steps with energy deposits (since non-depositing steps are
// 	skipped anyway).
//
//	trackID particleName X Y Z Ekin PX PY PZ
//	stepID Xi Yi Zi Ti KEi PXi PYi PZi Xf Yf Zf Tf KEf PXf PYf PZf
//
// <combiningLength> <unit>: Argument to /g4cmp/combiningLength, to have
//	StepAccumulator construct energy-weighted average position and
//	total energy of nearby steps.
//
// OR: testSecondaryProd -h to print this usage information
//
// 20220226  Michael Kelsey (TAMU)

#include "globals.hh"
#include "G4CMPSecondaryProduction.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPUtils.hh"
#include "G4Delete.hh"
#include "G4DynamicParticle.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UnitsTable.hh"
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <vector>


// Print help information

void Usage() {
  G4cout << "testSecondaryProd <trackFile> [combingingLength unit]" << G4endl
	 << "\nExercises G4CMPSecondaryProduction and G4CMPStepAccumulator"
	 << "\non a preset track with hits." << G4endl
	 << "<trackFile>: Name of text file describing all the steps of a"
	 << "\n   track; the first line describes the track kinematics, while"
	 << "\n   subesequent lines describe steps with energy deposits"
	 << "\n   (since non-depositing steps are skipped anyway)." << G4endl
	 << "\ntrackID particleName X Y Z Ekin PX PY PZ"
	 << "\nstepID Xi Yi Zi Ti KEi PXi PYi PZi Xf Yf Zf Tf KEf PXf PYf PZf"
	 << G4endl
	 << "\n<combiningLength> <unit>: Argument to /g4cmp/combiningLength,"
	 << "\n   to have StepAccumulator construct energy-weighted average"
	 << "\n   position and total energy of nearby steps."
	 << G4endl;
}


// Global variables for use in testing functions
// NOTE: It would be better to get rid of these, and pass arguments

namespace {
  G4CMPConfigManager* g4cmp = G4CMPConfigManager::Instance();
  G4LatticePhysical* lattice = 0;
  G4PVPlacement* volume = 0;
  G4CMPSecondaryProduction theProcess;
}

// Extract command-line arguments

G4bool GetArguments(int argc, char* argv[], G4String& trackFile,
		    G4double& combLen, G4String& units) {
  trackFile = units = ""; combLen = 0.;
  if (argc!=2 && argc!=4) return false;	// trackFile requires, others optional

  trackFile = argv[1];		// Must specify track definition file
  if (argc>2) {
    units = argv[3];
    combLen = strtod(argv[2], NULL) * G4UnitDefinition::GetValueOf(units);
  }

  return true;
}


// Single-volume geometry for tracking support

void ConstructGeometry() {
  G4String lname = "Ge";
  G4String mname = "G4_"+lname;

  // MUST USE 'new', SO THAT Geant4 CAN DELETE
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4Tubs* crystal = new G4Tubs("GeCrystal", 0., 5.*cm, 1.*cm, 0., 360.*deg);
  G4LogicalVolume* lv = new G4LogicalVolume(crystal, mat, crystal->GetName());
  volume = new G4PVPlacement(0,G4ThreeVector(), lv, lv->GetName(), 0,false,1);
  // NOTE: This PV is acting as the world volume
  // We will also need to create a G4TouchableHistory for it 

  lattice = G4LatticeManager::Instance()->LoadLattice(volume,lname);
}


// Read track information file and initialize track and steps list

G4bool InitializeTrack(std::istream& trackData, G4Track& theTrack) {
  if (!trackData.good()) return false;

  // Read track data from current (first) line of file

  return true;
}

G4bool LoadNextStep(std::istream& trackData,
		    std::map<G4int, G4Step>& theSteps) {
  if (!trackData.good()) return false;

  // NOTE: If input file ends with carriage return, may read a blank line

  return true;
}

G4bool LoadTrackInfo(const G4String& trackFile, G4Track& theTrack,
		     std::map<G4int, G4Step>& theSteps) {
  std::ifstream trackData(trackFile);
  if (!trackData) {
    G4cerr << "ERROR opening " << trackFile << G4endl; //TrackData to TraceFile?
    return false;
  }

  // Create track from first line of file, populating "vertex" data members
  if (!InitializeTrack(trackData, theTrack)) {
    G4cerr << "ERROR reading track data from " << trackFile << G4endl;
    return false;
  }

  // Populate vector of steps (not necessarily contiguous)
  while (trackData.good() && !trackData.eof()) {
    if (!LoadNextStep(trackData, theSteps)) {
      G4cerr << "ERROR reading step " << theSteps.size()+1 << " data from "
	     << trackFile << G4endl;
      return false;
    }
  }

  return true;
}


// Update track contents with information from step

void PrepareTrack(G4Track& theTrack, const std::pair<G4int,G4Step>& stepData) {
  G4int stepID = stepData.first;
  const G4Step& aStep = stepData.second;

  G4cout << " step " << theTrack.GetTrackID() << "/" << stepID
	 << " @ " << aStep.GetPostStepPoint()->GetPosition() << G4endl;
  
  // Copy post-step information into track

}


// Compare accumulated hit from process with actual step info

G4bool AnalyzeResults(const G4Track& theTrack, const G4Step& aStep) {
  G4bool resultGood = true;

  const G4CMPStepAccumulator* theHit = theProcess.GetAccumulator();

  // Verify that accumulator includes current step, and show content
  if (theHit->trackID == theTrack.GetTrackID() &&
      theHit->stepID == theTrack.GetCurrentStepNumber()) {
    G4cout << *theHit << G4endl;
  } else {
    G4cerr << "ERROR StepAccumulator does not match track state!" << G4endl
	   << "Accumulator step " << theHit->trackID << "/" << theHit->stepID
	   << G4endl;

    resultGood = false;
  }

  // Retrieve partition info, see if matches current step


  return resultGood;
}


// MAIN PROGRAM

int main(int argc, char* argv[]) {
  G4String trackFile, units;
  G4double combLen;
  if (!GetArguments(argc, argv, trackFile, combLen, units)) {
    G4cerr << "ERROR reading command line." << G4endl;
    Usage();
    return 1;
  }

  if (trackFile == "-h") {
    Usage();
    return 0;
  }

  G4cout << argv[0] << " processing " << trackFile;
  if (combLen > 0.) {
    G4cout << " with hit-combining " << G4BestUnit(combLen, "Length")
	   << G4endl;
  }

  // Load data from input file first
  G4Track theTrack;
  std::map<G4int, G4Step> theSteps;
  if (!LoadTrackInfo(trackFile, theTrack, theSteps)) return 2;

  // Set up the G4CMP framework, with tracking volume and lattice
  ConstructGeometry();

  theProcess.SetCombiningStepLength(combLen);

  G4cout << " processing track " << theTrack.GetTrackID() << " with "
	 << theSteps.size() << " recorded steps" << G4endl;

  G4VParticleChange* pchange=0;

  G4int nerrors = 0;
  for (const auto& stepData: theSteps) {
    PrepareTrack(theTrack, stepData);
    const G4Step& aStep = stepData.second;
    pchange = theProcess.PostStepDoIt(theTrack, aStep);
    if (!AnalyzeResults(theTrack, aStep)) nerrors++;
  }

  G4cout << "Processing complete with " << nerrors << " errors" << G4endl;

  return nerrors;
}
