////////////////////////////////////////////////////////////////////
//
// ValidationAnalysis.cc
// linehan3@fnal.gov
//
// This script takes stepping output from running G4CMP
// and organizes it for easy analysis. 
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"


//---------------------------------------------------------------------------
// Define a set of structs for use interpreting the output from G4CMP
struct Step {
  int runID;
  int eventID;
  int trackID;
  std::string particleName;
  double preStepEnergy_eV;
  double preStepKinEnergy_eV;
  double preStepX_mm;
  double preStepY_mm;
  double preStepZ_mm;
  double preStepT_ns;

  double postStepEnergy_eV;
  double postStepKinEnergy_eV;
  double postStepX_mm;
  double postStepY_mm;
  double postStepZ_mm;
  double postStepT_ns;
  size_t nReflections;
  std::string processName;
};

struct Track {
  int runID;
  int eventID;
  int trackID;
  std::vector<Step> stepVect;
};

struct Event {
  int runID;
  int eventID;
  std::vector<Track> trackVect;
};

//----------------------------------------------------------------------------
// Parsing function
std::vector<Event> ReadInG4CMPStepTextFile(std::string filename) {
  std::vector<Event> output;
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  int eventID = -1;
  int runID = -1;
  int trackID = -1;
  int counter = 0;
  while (1) {
    if (!infile.good()) break;
    if (infile.is_open()) {
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between commas)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while (getline(check1,token,' ')) {
	tokens.push_back(token);
      }
      if (tokens.size() == 0) break;


      //If our event or run number changes, then start a new event and push it
      //back into the output
      bool startNewTrack = false;
      if (std::atoi(tokens[0].c_str()) != runID || std::atoi(tokens[1].c_str()) != eventID) {
	Event newEvent;
	newEvent.runID = std::atoi(tokens[0].c_str());
	newEvent.eventID = std::atoi(tokens[1].c_str());
	output.push_back(newEvent);
	runID = newEvent.runID;
	eventID = newEvent.eventID;
	startNewTrack = true;
	
	if (eventID % 10 == 0) std::cout << "Done reading " << eventID << " events." << std::endl;
	counter++;
      }

      
      //If our track number changes (or if we roll over to a new event), then start a new track and push back into the current event
      if (std::atoi(tokens[2].c_str()) != trackID || startNewTrack) {
	Track newTrack;
	newTrack.trackID = std::atoi(tokens[2].c_str());
	newTrack.eventID = eventID;
	newTrack.runID = runID;
	trackID = newTrack.trackID;
	output[output.size()-1].trackVect.push_back(newTrack);
      }

      //Log the step information
      Step newStep;
      newStep.runID = std::atoi(tokens[0].c_str());
      newStep.eventID = std::atoi(tokens[1].c_str());
      newStep.trackID = std::atoi(tokens[2].c_str());
      newStep.particleName = tokens[3];
      
      newStep.preStepX_mm = std::atof(tokens[4].c_str());
      newStep.preStepY_mm = std::atof(tokens[5].c_str());
      newStep.preStepZ_mm = std::atof(tokens[6].c_str());
      newStep.preStepT_ns = std::atof(tokens[7].c_str());      
      newStep.preStepEnergy_eV = std::atof(tokens[8].c_str());
      newStep.preStepKinEnergy_eV = std::atof(tokens[9].c_str());
      
      newStep.postStepX_mm = std::atof(tokens[10].c_str());
      newStep.postStepY_mm = std::atof(tokens[11].c_str());
      newStep.postStepZ_mm = std::atof(tokens[12].c_str());
      newStep.postStepT_ns = std::atof(tokens[13].c_str());
      newStep.postStepEnergy_eV = std::atof(tokens[14].c_str());
      newStep.postStepKinEnergy_eV = std::atof(tokens[15].c_str());
      
      newStep.nReflections = std::atoi(tokens[16].c_str());
      
      newStep.processName = tokens[17];
      
      int currentTLL = output[output.size()-1].trackVect.size();
      output[output.size()-1].trackVect[currentTLL-1].stepVect.push_back(newStep);
    }
  }
  return output;
}


//----------------------------------------------------------------------------
// Analysis of boundary transmission validation test
// (input = Validation_BoundaryTransmission.txt)
void ValidationAnalysis_BoundaryTransmission() {
  //Read in the events from the textfile
  std::vector<Event> eventVect
    = ReadInG4CMPStepTextFile("Validation_BoundaryTransmission.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge
  //interface
  TFile * outFile
    = new TFile("Validation_BoundaryTransmission.root","RECREATE");

  //Histogram definitions
  double diskRadius_mm = 10.0;
  double z_SiGeInterface_mm = -0.5;
  double maxZ_mm = 1.1;
  double minZ_mm = -2.1;
  TH1F * h_boundaryTransmission =
    new TH1F("h_boundaryTransmission",
	     "h_boundaryTransmission; All/Reflected/Transmitted; Surface-Incident Steps",3,0,3);
  
  TH2F * h_postStepPointXY =
    new TH2F("h_postStepPointXY","h_postStepPointXY; Step X [mm]; Step Y [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);
  TH2F * h_postStepPointXZ =
    new TH2F("h_postStepPointXZ","h_postStepPointXZ; Step X [mm]; Step Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  TH2F * h_postStepPointYZ =
    new TH2F("h_postStepPointYZ","h_postStepPointYZ; Step Y [mm]; Step Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);				      				      
  TH2F * h_reflectedFromAboveNextStepEndingsXY =
    new TH2F("h_reflectedFromAboveNextStepEndingsXY",
	     "h_reflectedFromAboveNextStepEndingsXY; Next step ending X [mm]; Next step ending Y [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);

  TH2F * h_reflectedFromAboveNextStepEndingsXZ =
    new TH2F("h_reflectedFromAboveNextStepEndingsXZ",
	     "h_reflectedFromAboveNextStepEndingsXZ; Next step ending X [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  TH2F * h_reflectedFromAboveNextStepEndingsYZ =
    new TH2F("h_reflectedFromAboveNextStepEndingsYZ",
	     "h_reflectedFromAboveNextStepEndingsYZ; Next step ending Y [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_transmittedFromAboveNextStepEndingsXY =
    new TH2F("h_transmittedFromAboveNextStepEndingsXY",
	     "h_transmittedFromAboveNextStepEndingsXY; Next step ending X [mm]; Next step ending Y [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);
  
  TH2F * h_transmittedFromAboveNextStepEndingsXZ =
    new TH2F("h_transmittedFromAboveNextStepEndingsXZ",
	     "h_transmittedFromAboveNextStepEndingsXZ; Next step ending X [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
 
  TH2F * h_transmittedFromAboveNextStepEndingsYZ =
    new TH2F("h_transmittedFromAboveNextStepEndingsYZ",
	     "h_transmittedFromAboveNextStepEndingsYZ; Next step ending Y [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_reflectedFromBelowNextStepEndingsXY =
    new TH2F("h_reflectedFromBelowNextStepEndingsXY",
	     "h_reflectedFromBelowNextStepEndingsXY; Next step ending X [mm]; Next step ending Y [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);
  
  TH2F * h_reflectedFromBelowNextStepEndingsXZ =
    new TH2F("h_reflectedFromBelowNextStepEndingsXZ",
	     "h_reflectedFromBelowNextStepEndingsXZ; Next step ending X [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_reflectedFromBelowNextStepEndingsYZ =
    new TH2F("h_reflectedFromBelowNextStepEndingsYZ",
	     "h_reflectedFromBelowNextStepEndingsYZ; Next step ending Y [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_transmittedFromBelowNextStepEndingsXY =
    new TH2F("h_transmittedFromBelowNextStepEndingsXY",
	     "h_transmittedFromBelowNextStepEndingsXY; Next step ending X [mm]; Next step ending Y [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);
  
  TH2F * h_transmittedFromBelowNextStepEndingsXZ =
    new TH2F("h_transmittedFromBelowNextStepEndingsXZ",
	     "h_transmittedFromBelowNextStepEndingsXZ; Next step ending X [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_transmittedFromBelowNextStepEndingsYZ =
    new TH2F("h_transmittedFromBelowNextStepEndingsYZ",
	     "h_transmittedFromBelowNextStepEndingsYZ; Next step ending Y [mm]; Next step ending Z [mm]",
	     1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over pairs of steps within tracks
      for (int iS = 0; iS < theTrack.stepVect.size()-2; ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
	Step theStep2 = theTrack.stepVect[iS+1];
	Step theStep3 = theTrack.stepVect[iS+2];
       	
	//Compute quantities
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	
	double x0_2 = theStep2.preStepX_mm;
	double y0_2 = theStep2.preStepY_mm;
	double z0_2 = theStep2.preStepZ_mm;
	double x1_2 = theStep2.postStepX_mm;
	double y1_2 = theStep2.postStepY_mm;
	double z1_2 = theStep2.postStepZ_mm;
	double dX_2 = x1_2-x0_2;
	double dY_2 = y1_2-y0_2;
	double dZ_2 = z1_2-z0_2;
	double r1_2 = TMath::Power(x1_2*x1_2 + y1_2*y1_2,0.5);
	double stepLength2_mm =
	  TMath::Power(dX_2*dX_2 + dY_2*dY_2 + dZ_2*dZ_2,0.5 );
	
	double dz_2 = z1_2-z0_2;
	double z_dir2 = dz_2 / fabs(dz_2);
	
	double x0_3 = theStep3.preStepX_mm;
	double y0_3 = theStep3.preStepY_mm;
	double z0_3 = theStep3.preStepZ_mm;
	double x1_3 = theStep3.postStepX_mm;
	double y1_3 = theStep3.postStepY_mm;
	double z1_3 = theStep3.postStepZ_mm;
	double dX_3 = x1_3-x0_3;
	double dY_3 = y1_3-y0_3;
	double dZ_3 = z1_3-z0_3;
	double r1_3 = TMath::Power(x1_3*x1_3 + y1_3*y1_3,0.5);
	double stepLength3_mm =
	  TMath::Power(dX_3*dX_3 + dY_3*dY_3 + dZ_3*dZ_3,0.5 );
	
	double dz_3 = z1_3-z0_3;
	double z_dir3 = dz_3 / fabs(dz_3);
	
	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);

	//First conditional, to get a pure set of internal-boundary-hitters:
	//1. If the first step is long and ends within 1 nm of the Si-Ge
	//interface and not on the walls, fill a histo with an entry (i.e.
	//all "incidents")
	//2. If the next step is short and the next-next step goes in the
	//opposite z-direction to the first, then this is a reflection. Fill
	//the same histo with a different entry representing reflection.
	//3. If the next step is long and the next-next step goes in the same
	//z-direction to the first, then this is a transmission. Fill the same
	//histo with a different entry representing transmission.
	if (stepLength1_mm > 0.001) {

	  //All boundary hitters
	  if (r1 < (0.99*diskRadius_mm) && fabs(z1-z_SiGeInterface_mm) < 1e-9) {
	    h_boundaryTransmission->Fill(0);

	    //Looking at next-step and next-next step. If directions are not
	    //the same, then it's a reflection
	    if (stepLength2_mm < 1E-10 && z_dir1 / z_dir3 == -1) {
	      h_boundaryTransmission->Fill(1);

	      //Some cross-checks
	      if (z0 > z_SiGeInterface_mm) {
		h_reflectedFromAboveNextStepEndingsXY->Fill(x1_3,y1_3);
		h_reflectedFromAboveNextStepEndingsXZ->Fill(x1_3,z1_3);
		h_reflectedFromAboveNextStepEndingsYZ->Fill(y1_3,z1_3);
	      } else {
		h_reflectedFromBelowNextStepEndingsXY->Fill(x1_3,y1_3);
		h_reflectedFromBelowNextStepEndingsXZ->Fill(x1_3,z1_3);
		h_reflectedFromBelowNextStepEndingsYZ->Fill(y1_3,z1_3);
	      }
	    } else {
	      //Otherwise, it's a transmission
	      h_boundaryTransmission->Fill(2);
	      
	      //Some cross-checks
	      if (z0 > z_SiGeInterface_mm) {
		h_transmittedFromAboveNextStepEndingsXY->Fill(x1_2,y1_2);
		h_transmittedFromAboveNextStepEndingsXZ->Fill(x1_2,z1_2);
		h_transmittedFromAboveNextStepEndingsYZ->Fill(y1_2,z1_2);
	      } else {
		h_transmittedFromBelowNextStepEndingsXY->Fill(x1_2,y1_2);
		h_transmittedFromBelowNextStepEndingsXZ->Fill(x1_2,z1_2);
		h_transmittedFromBelowNextStepEndingsYZ->Fill(y1_2,z1_2);
	      }
	    }
	  }
	}
      }
    }
  }
  outFile->Write();
}



//--------------------------------------------------------------------------
// Analysis of polycrystalline elastic scattering validation test
// (input = Validation_PolycrystalElasticScattering.txt)
void ValidationAnalysis_PolycrystalElasticScattering()
{
  //Read in the events from the textfile
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile("Validation_PolycrystalElasticScattering.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge interface
  TFile * outFile = new TFile("Validation_PolycrystalElasticScattering.root","RECREATE");

  
  //Histogram definitions
  double diskRadius_mm = 1e-3;
  double minZ_mm = -1.500110;
  double maxZ_mm = -1.499990;
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","h_postStepPointXY; Step X [mm]; Step Y [mm]",
				      1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
				      1000,-diskRadius_mm*1.1,diskRadius_mm*1.1);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","h_postStepPointXZ; Step X [mm]; Step Z [mm]",
				      1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
				      1000,minZ_mm,maxZ_mm);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","h_postStepPointYZ; Step Y [mm]; Step Z [mm]",
				      1000,-diskRadius_mm*1.1,diskRadius_mm*1.1,
				      1000,minZ_mm,maxZ_mm);
  TH1F * h_stepLength = new TH1F("h_stepLength","h_stepLength; Step Length [mm]; Steps",1000,0,0.001);
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
       	
	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);
	
	//Step length
	h_stepLength->Fill(stepLength1_mm);
      }
    }
  }
  outFile->Write();
}

//-----------------------------------------------------------------------------
// Analysis of SC pairbreaking validation test
// (input = Validation_SCPairbreaking.txt)
void ValidationAnalysis_SCPairbreaking() {

  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_SCPairbreaking.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge
  //interface
  TFile * outFile = new TFile("Validation_SCPairbreaking.root","RECREATE");
  
  //Histogram definitions. Note the boilerplate postStepPoint plots are only
  //for phonons
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  TH2F * h_postStepPointXY =
    new TH2F("h_postStepPointXY","h_postStepPointXY; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ =
    new TH2F("h_postStepPointXZ","h_postStepPointXZ; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ =
    new TH2F("h_postStepPointYZ","h_postStepPointYZ; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_stepDurationVsPhononEnergy =
    new TH2F("h_stepDurationVsPhononEnergy",
	     "h_stepDurationVsPhononEnergy; Pre-Step Energy [meV]; Step Duration [ns]",
	     nStepsE,minE_meV,maxE_meV,
	     1000,0,1.0);

  TH2F * h_qpEnergyVsPreStepEnergy =
    new TH2F("h_qpEnergyVsPreStepEnergy",
	     "QP Energy vs. Pre-Step (Phonon) Energy; Phonon energy [meV]; QP Energies [meV]",
	     nStepsE*10,minE_meV,maxE_meV,
	     nStepsE*10,minE_meV,maxE_meV);
  
  //Loop over events
  double thisPhononEnergy_eV = -1;
  int thisPhononEvent = -1;
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	
	//Select only phonons
	if (particleName.find("phonon") != std::string::npos) {

	  //Fill cross-check histograms
	  h_postStepPointXY->Fill(x1,y1);
	  h_postStepPointXZ->Fill(x1,z1);
	  h_postStepPointYZ->Fill(y1,z1);
	  
	  //Now select only phonons which have scPairbreaking as the process
	  if (processName.find("scPairBreaking") != std::string::npos) {
	    
	    //Plot the step duration vs. the phonon energy
	    h_stepDurationVsPhononEnergy->Fill(preStepEnergy_eV*1000,t1);

	    //Also log the QP energy
	    thisPhononEnergy_eV = preStepEnergy_eV;
	    thisPhononEvent = theStep1.eventID;
	  }
	}

	//Next, select QPs from the same event
	if (particleName.find("QP") != std::string::npos) {

	  //Select only first steps of the QP
	  if (iS == 0) {
	    
	    //Check to make sure that this QP's event is the same as the phonon
	    //we just extracted an energy from	  
	    if (theStep1.eventID == thisPhononEvent) {
	      
	      //Fill cross-check histograms. Have to use kinetic energy because
	      //the QPs are given a mass energy that is artificial so that they
	      //get a velocity that is acceptable for diffusion physics.
	      h_qpEnergyVsPreStepEnergy->Fill(thisPhononEnergy_eV*1000.,
					      theStep1.preStepKinEnergy_eV*1000.);
	    }
	  }
	}
      }
    }
  }

  //Now do fitting of each of the slices of the time histogram and make a graph
  //out of it. Setup for fitter and other misc things
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000);  
  double fitMin_ns = 0.01;
  //double fitMax_ns = 0.4; //for low energy
  double fitMax_ns = 0.2; //For high energy
  double gap_meV = 0.176;
  double tau0_ph_ns = 0.242;
   
  //Establish the graph of lifetimes
  TGraphErrors * g_fitExponents = new TGraphErrors();

  //Do the actual fits
  int offset = 0;
  TF1 * expFit =
    new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",fitMin_ns,fitMax_ns);
  
  for (int iBX = 1; iBX < h_stepDurationVsPhononEnergy->GetNbinsX(); ++iBX) {

    //Check to see if we are below 2Delta. If so, skip
    if (h_stepDurationVsPhononEnergy->GetXaxis()->GetBinCenter(iBX) <
	2*gap_meV) {
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }

    //Otherwise, don't.
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice =
      (TH1F*)h_stepDurationVsPhononEnergy->ProjectionY(name,iBX,iBX,"");
    
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");
    std::cout << "iBX: " << iBX << std::endl;

    if (iBX < ((int)h_stepDurationVsPhononEnergy->GetNbinsX() / 2.0 )) {

      //Hardcoded start point for low energy 
      expFit->SetParameter(0,2000./(theSlice->GetMean()*60));

      //Hardcoded start point for low-energy
      expFit->SetParameter(1,0.19); 
    } else if (iBX == 20) {      
      expFit->SetParameter(0,200); //Hardcoded start point (based on # entries)
      expFit->SetParameter(1,0.15); //Hardcoded start point for high-energy
    } else {
      expFit->SetParameter(0,2000); //Hardcoded start point (based on # entries)
      expFit->SetParameter(1,0.11); //Hardcoded start point for high-energy
    }
    theSlice->Fit(expFit,"","",fitMin_ns,fitMax_ns);
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);
    std::cout << "fitExponent_ns: " << fitExponent_ns << std::endl;
    
    //If the fit exponent value is less than some value, then retry    
    //Convert the fit exponent to tau0/tauB
    double tau0DivTauB = tau0_ph_ns / fitExponent_ns;
    double tau0DivTauB_error = tau0DivTauB *
      (fitExponentError_ns / fitExponent_ns);    
    g_fitExponents->SetPoint(iBX-offset-1,h_stepDurationVsPhononEnergy->GetXaxis()->GetBinCenter(iBX),tau0DivTauB);
    g_fitExponents->SetPointError(iBX-offset-1,0,tau0DivTauB_error);
  }
  g_fitExponents->SetName("g_fitExponents");
  g_fitExponents->SetTitle("Fit Exponents");
  g_fitExponents->GetXaxis()->SetTitle("Pre-Step Energy [meV]");
  g_fitExponents->GetYaxis()->SetTitle("Tau0 / TauB");
  g_fitExponents->Draw("ALP");  
  g_fitExponents->Write();
  outFile->Write();
}


//----------------------------------------------------------------------------
// Analysis of QP Phonon Radiation
// (input = Validation_QPPhononRadiation.txt)
void ValidationAnalysis_QPPhononRadiation() {
  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_QPPhononRadiation.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute phonon radiation properies
  TFile * outFile = new TFile("Validation_QPPhononRadiation.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointXY_qpRad =
    new TH2F("h_postStepPointXY_qpRad",
	     "h_postStepPointXY, Radiation Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpRad =
    new TH2F("h_postStepPointXZ_qpRad",
	     "h_postStepPointXZ, Radiation Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpRad =
    new TH2F("h_postStepPointYZ_qpRad",
	     "h_postStepPointYZ, Radiation Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_radiatedEnergyVsQPPreStepEnergy =
    new TH2F("h_radiatedEnergyVsQPPreStepEnergy",
	     "h_radiatedEnergyVsQPPreStepEnergy; QP Pre-Step Energy [meV]; Energy Radiated [meV];",
	     10*nStepsE,minE_meV,maxE_meV,
	     10*nStepsE,minE_meV,maxE_meV);

  TH2F * h_remainingEnergyVsQPPreStepEnergy =
    new TH2F("h_remainingEnergyVsQPPreStepEnergy",
	     "h_remainingEnergyVsQPPreStepEnergy; QP Pre-Step Energy [meV]; Energy Radiated [meV];",
	     10*nStepsE,minE_meV,maxE_meV,
	     10*nStepsE,minE_meV,maxE_meV);
  
  TH2F * h_stepDurationVsQPPreStepEnergy =
    new TH2F("h_stepDurationVsQPPreStepEnergy",
	     "h_stepDurationVsQPPreStepEnergy; QP Pre-Step Energy [meV]; Time Since Last Radiation [ns]",
	     nStepsE,minE_meV,maxE_meV*0.22,
	     10000,0,1000);
    
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);
	  
	  //Now select only phonons which have scPairbreaking as the process
	  if (processName.find("qpRadiatesPhonon") != std::string::npos) {

	    //For this track, compute the difference in time between the last
	    //radiation step and the current radiation step,
	    //and then make this the last radiation step.
	    double thisDeltaT_ns = t1-lastRadStepTime_ns;
	    lastRadStepTime_ns = t1;

	    //Plot positions (now physically meaningful)
	    h_postStepPointXY_qpRad->Fill(x1,y1);
	    h_postStepPointXZ_qpRad->Fill(x1,z1);
	    h_postStepPointYZ_qpRad->Fill(y1,z1);

	    //Plot the deltaT since the last radiation step against the QP
	    //pre-energy
	    h_stepDurationVsQPPreStepEnergy->Fill(preStepKinEnergy_eV*1000,
						  thisDeltaT_ns);

	    
	    //Plot the QP's deltaE (i.e. the phonon's radiated energy) vs. QP
	    //pre-energy
	    double deltaE_eV = preStepKinEnergy_eV - postStepKinEnergy_eV;
	    h_radiatedEnergyVsQPPreStepEnergy->Fill(preStepKinEnergy_eV*1000,
						    deltaE_eV*1000);
	    h_remainingEnergyVsQPPreStepEnergy->Fill(preStepKinEnergy_eV*1000,
						     postStepKinEnergy_eV*1000);
	  }
	}
      }
    }
  }
  
  //Now do fitting of each of the slices of the time histogram and make a graph
  //out of it
  //Set up fitter to do lots of iterations
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000); 
  double fitMin_ns = 0.1;
  double fitMax_ns = 500; //For high energy
  double gap_meV = 0.176;
  double tau0_qp_ns = 438.;
  
  //Establish the graph of lifetimes
  TGraphErrors * g_fitExponents = new TGraphErrors();

  //Now do the fits
  int offset = 0;
  TF1 * expFit =
    new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",fitMin_ns,fitMax_ns);
  
  for (int iBX = 1; iBX < h_stepDurationVsQPPreStepEnergy->GetNbinsX(); ++iBX) {
    std::cout << " pre-step energy bin center: "
	      << h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX)
	      << std::endl;
    if (h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX)
	< 1*gap_meV) {
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }
    
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice =
      (TH1F*)h_stepDurationVsQPPreStepEnergy->ProjectionY(name,iBX,iBX,"");
    
    //if( h_stepDurationVsQPPreStepEnergy->GetBinCenter(iBX)
    //< 4*gap_meV ){ theSlice->Rebin(5); }
    
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");

    //Hardcoded start point for low energy 
    //expFit->SetParameter(0,2000./(theSlice->GetMean()*60)); 
    //expFit->SetParameter(1,0.35); //Hardcoded start point for low-energy

    //Hardcoded start point (based on number of entries)
    expFit->SetParameter(0,theSlice->GetBinContent(1));
    //Hardcoded start point for high-energy
    expFit->SetParameter(1,theSlice->GetMean());
    theSlice->Fit(expFit,"","",fitMin_ns,fitMax_ns);
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);

    //Convert the fit exponent to taus/tau0
    double tauSDivTau0 = fitExponent_ns / tau0_qp_ns;
    double tauSDivTau0_error = tauSDivTau0 *
      (fitExponentError_ns / fitExponent_ns);
    g_fitExponents->SetPoint(iBX-offset-1,h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX),tauSDivTau0);
    g_fitExponents->SetPointError(iBX-offset-1,0,tauSDivTau0_error);
  }
  g_fitExponents->SetName("g_fitExponents");
  g_fitExponents->SetTitle("Fit Exponents");
  g_fitExponents->GetXaxis()->SetTitle("Pre-Step Energy [meV]");
  g_fitExponents->GetYaxis()->SetTitle("TauS / Tau0");
  
  TCanvas * c0 = new TCanvas();
  g_fitExponents->Draw("ALP");  
  g_fitExponents->Write();
  outFile->Write();
}




//------------------------------------------------------------------------------
// Analysis of QP Recombination
// (input = Validation_QPRecombination.txt)
void ValidationAnalysis_QPRecombination() {
  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_QPRecombination.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute QP recombination properties
  TFile * outFile = new TFile("Validation_QPRecombination.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_qpRecomb =
    new TH2F("h_postStepPointXY_qpRecomb",
	     "h_postStepPointXY, Recombination Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpRecomb =
    new TH2F("h_postStepPointXZ_qpRecomb",
	     "h_postStepPointXZ, Recombination Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpRecomb =
    new TH2F("h_postStepPointYZ_qpRecomb",
	     "h_postStepPointYZ, Recombination Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_stepDurationVsQPPreStepEnergy =
    new TH2F("h_stepDurationVsQPPreStepEnergy",
	     "h_stepDurationVsQPPreStepEnergy; QP Pre-Step Energy [meV]; Time Since QP Creation [ns]",
	     nStepsE,minE_meV,maxE_meV*0.22,
	     50000,0,5e6); //For 0.2K
  //50000,0,5e3); //For 0.5K

  TH2F * h_phononEnergyVsQPEnergy =
    new TH2F("h_phononEnergyVsQPEnergy",
	     "h_phononEnergyVsQPEnergy; QP Pre-Step Energy [meV]; Resultant Phonon Energy [meV]",
	     nStepsE*10,minE_meV,maxE_meV,
	     nStepsE*10,minE_meV,maxE_meV);
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    double thisEventQPEnergy = -1;
    int thisEvent = -1;
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);
	  
	  //Now select only phonons which have scPairbreaking as the process
	  if (processName.find("qpRecombination") != std::string::npos) {

	    //Plot positions (now physically meaningful)
	    h_postStepPointXY_qpRecomb->Fill(x1,y1);
	    h_postStepPointXZ_qpRecomb->Fill(x1,z1);
	    h_postStepPointYZ_qpRecomb->Fill(y1,z1);

	    //Plot the deltaT since the last radiation step against the QP
	    //pre-energy
	    h_stepDurationVsQPPreStepEnergy->Fill(preStepKinEnergy_eV*1000,t1);

	    //Note the QP's energy at recombination to compare to phonon
	    thisEventQPEnergy = preStepKinEnergy_eV;
	    thisEvent = theStep1.eventID;
	  }	  
	}

	//Select only phonons
	if (particleName.find("phonon") != std::string::npos) {

	  //Make sure we're looking at the same event
	  if (theStep1.eventID == thisEvent) {

	    //Only fill once for each phonon
	    if (iS == 0) {
	      h_phononEnergyVsQPEnergy->Fill(thisEventQPEnergy*1000,preStepEnergy_eV*1000);
	    }
	  }
	}
      }
    }
  }
  
  //Now do fitting of each of the slices of the time histogram and make a
  //graph out of it
  //Set up fitter to do lots of iterations
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000); 
  double fitMin_ns = 0.1;
  //  double fitMax_ns = 300e5; //For high energy
  double gap_meV = 0.176;
  double tau0_qp_ns = 438.;
  
  //Establish the graph of lifetimes
  TGraphErrors * g_fitExponents = new TGraphErrors();

  //Now do the fits
  int offset = 0;
  for (int iBX = 1; iBX < h_stepDurationVsQPPreStepEnergy->GetNbinsX(); ++iBX) {
    std::cout << " pre-step energy bin center: "
	      << h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX)
	      << std::endl;
    if (h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX) <
	1*gap_meV) {
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }
    
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice =
      (TH1F*)h_stepDurationVsQPPreStepEnergy->ProjectionY(name,iBX,iBX,"");

    //Want to put expfit here because that way we can adaptively select the
    //slice mean
    TF1 * expFit = new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",
			   fitMin_ns,2*theSlice->GetMean());
    
    //if( h_stepDurationVsQPPreStepEnergy->GetBinCenter(iBX)
    //< 4*gap_meV ){ theSlice->Rebin(5); }
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");
    theSlice->Rebin(20);

    //Hardcoded start point for low energy 
    //expFit->SetParameter(0,2000./(theSlice->GetMean()*60));
    //expFit->SetParameter(1,0.35); //Hardcoded start point for low-energy

    //Hardcoded start point (based on number of entries)
    expFit->SetParameter(0,theSlice->GetBinContent(1));
    //Hardcoded start point for high-energy
    expFit->SetParameter(1,theSlice->GetMean());
    //expFit->SetParameter(1,0.5e6); //Hardcoded start point for high-energy
    theSlice->Fit(expFit,"","",fitMin_ns,2*theSlice->GetMean());
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);

    //Convert the fit exponent to taur/tau0
    double tauRDivTau0 = fitExponent_ns / tau0_qp_ns;
    double tauRDivTau0_error =
      tauRDivTau0 * (fitExponentError_ns / fitExponent_ns);
    g_fitExponents->SetPoint(iBX-offset-1,h_stepDurationVsQPPreStepEnergy->GetXaxis()->GetBinCenter(iBX),tauRDivTau0);
    g_fitExponents->SetPointError(iBX-offset-1,0,tauRDivTau0_error);
  }
  g_fitExponents->SetName("g_fitExponents");
  g_fitExponents->SetTitle("Fit Exponents");
  g_fitExponents->GetXaxis()->SetTitle("Pre-Step Energy [meV]");
  g_fitExponents->GetYaxis()->SetTitle("TauR / Tau0");
  
  TCanvas * c0 = new TCanvas();
  g_fitExponents->Draw("ALP");  
  g_fitExponents->Write();
  outFile->Write();
}

//-----------------------------------------------------------------------------
// Calculates the time of first passage density function for the interval.
// This is bidirectional (assumes QP-killing BC on both ends). This function
// is used  in the validation function below
TGraph * CalculateTOFPOnTheInterval() {
  int nPts = 10000;
  TGraph * g_lowTLimit = new TGraph(nPts);
  TGraph * g_highTLimit = new TGraph(nPts);
  TGraph * g_out = new TGraph(nPts);
  
  double min_t = 0.001;
  double max_t = 2;
  double delta_t = (max_t - min_t) / ((double)nPts);
  double sumForNormalization = 0;
  for (int iP = 0; iP < nPts; ++iP) {
    double t = delta_t*iP + min_t;
    double lowEVal
      = 8/TMath::Sqrt(TMath::Pi()) *
      (TMath::Exp(-1.0/16.0/t) / 16 / t / t *
       (TMath::Sqrt(t) - 8*TMath::Power(t,1.5) + 192*TMath::Power(t,2.5) ) +
       TMath::Exp(-1.0/16.0/t)*
       (0.5*TMath::Power(t,-0.5) - 12*TMath::Power(t,0.5)
	+ 480*TMath::Power(t,1.5)));
    
    double highEVal = 4*TMath::Pi()*TMath::Exp(-TMath::Pi()*TMath::Pi()*t);
    g_lowTLimit->SetPoint(iP,t,lowEVal);
    g_highTLimit->SetPoint(iP,t,highEVal);
    
    
    //Match limits
    //Scaling to get us to real time (instead of dimensionless time).
    double t_scaling_ns = 2.0*2.0/18.0*1e9; 
    // Above, 2 cm is used as it is the length of the bar.
    // Above, 18 cm2/s is used as it is the diffusion constant expected for
    //near-gap QPs, given energy dependence.
    if (t < 0.033) {      
      g_out->SetPoint(iP,t*t_scaling_ns,lowEVal);
      sumForNormalization += (lowEVal * delta_t * t_scaling_ns);
    } else {
      g_out->SetPoint(iP,t*t_scaling_ns,highEVal);
      sumForNormalization += (highEVal * delta_t * t_scaling_ns);
    }
  }
  
  TCanvas * c1 = new TCanvas();
  g_lowTLimit->SetLineColor(kBlue);
  g_highTLimit->SetLineColor(kRed);
  g_lowTLimit->Draw("AL");
  g_highTLimit->Draw("Lsame");
  g_lowTLimit->GetXaxis()->SetTitle("t [dimensionless]");
  g_lowTLimit->GetYaxis()->SetTitle("#rho(t)");
  //  sumForNormalization = 450;
  g_out->Scale(1.0/((double)sumForNormalization));
  g_out->SetLineColor(kBlack);
  g_out->SetLineStyle(9);
  g_out->Draw("Lsame");
  
  TLegend * l1 = new TLegend();
  l1->AddEntry(g_lowTLimit,"Low-T limit");
  l1->AddEntry(g_highTLimit,"High-T limit");
  l1->AddEntry(g_out,"Combined");
  l1->Draw("same");
  
  return g_out;  
}


//-----------------------------------------------------------------------------
// Analysis of QP Diffusion (Time-of-first-passage checks)
// (input = Validation_QPDiffusionTOFP)
void ValidationAnalysis_QPDiffusionTOFP() {
  //Read in the events from the textfile
  std::vector<Event> eventVect
    = ReadInG4CMPStepTextFile("Validation_QPDiffusionTOFP.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the
  //Si-Ge interface
  TFile * outFile = new TFile("Validation_QPDiffusionTOFP.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  double minTime_ns = 0;
  double maxTime_ns = 3e8;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointXY_qpLast =
    new TH2F("h_postStepPointXY_qpLast",
	     "h_postStepPointXY, Last QP Steps; Step X [mm]; Step Y [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpLast =
    new TH2F("h_postStepPointXZ_qpLast",
	     "h_postStepPointXZ, Last QP Steps; Step X [mm]; Step Z [mm]",
	     1000,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpLast =
    new TH2F("h_postStepPointYZ_qpLast",
	     "h_postStepPointYZ, Last QP Steps; Step Y [mm]; Step Z [mm]",
	     1000,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH1F * h_timeOfLastStep =
    new TH1F("h_timeOfLastStep","h_timeOfLastStep; Time [ns]; QPs;",
	     10000,minTime_ns,maxTime_ns);

  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);

	  //Only fill this if the step is the last one
	  if (iS == theTrack.stepVect.size()-1) {
	    h_postStepPointXY_qpLast->Fill(x1,y1);
	    h_postStepPointXZ_qpLast->Fill(x1,z1);
	    h_postStepPointYZ_qpLast->Fill(y1,z1);
	    
	    h_timeOfLastStep->Fill(t1);
	  }	  
	}
      }
    }
  }

  //Let's compare to a theory curve
  TGraph * g_TOFP_expected = CalculateTOFPOnTheInterval();  
  g_TOFP_expected->Write();

  TCanvas * c10 = new TCanvas();
  h_timeOfLastStep->Draw("HIST");
  g_TOFP_expected->SetLineStyle(1);
  g_TOFP_expected->SetLineColor(kRed);
  double binWidth_ns = h_timeOfLastStep->GetBinWidth(1);
  h_timeOfLastStep->Scale(1.0/((double)h_timeOfLastStep->Integral())/binWidth_ns);
  h_timeOfLastStep->Draw("HIST");
  g_TOFP_expected->Draw("Lsame");
  TLegend * l2 = new TLegend();
  l2->AddEntry(h_timeOfLastStep,"TOFP, Either Nb End");
  l2->AddEntry(g_TOFP_expected,
	       "TOFP, Theoretically Expected for Sim. Geometry");
  
  l2->Draw("same");
  c10->SaveAs("TOFPComparison.png");
  
  outFile->Write();
}


//----------------------------------------------------------------------------
// Analysis of QP Diffusion (Recombination death locations)
// (input = Validation_QPDiffusionPlusRecombination)
// This is literally the same analysis as the previous one, minus the TOFP bit.
void ValidationAnalysis_QPDiffusionPlusRecombination() {
  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_QPDiffusionPlusRecombination.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge
  //interface
  TFile * outFile =
    new TFile("Validation_QPDiffusionPlusRecombination.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  double minTime_ns = 0;
  double maxTime_ns = 3e8;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_qpLast =
    new TH2F("h_postStepPointXY_qpLast",
	     "h_postStepPointXY, Last QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpLast =
    new TH2F("h_postStepPointXZ_qpLast",
	     "h_postStepPointXZ, Last QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpLast =
    new TH2F("h_postStepPointYZ_qpLast",
	     "h_postStepPointYZ, Last QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH1F * h_timeOfLastStep = new TH1F("h_timeOfLastStep",
				     "h_timeOfLastStep; Time [ns]; QPs;",
				     10000,minTime_ns,maxTime_ns);
  
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);

	  //Only fill this if the step is the last one
	  if (iS == theTrack.stepVect.size()-1) {
	    h_postStepPointXY_qpLast->Fill(x1,y1);
	    h_postStepPointXZ_qpLast->Fill(x1,z1);
	    h_postStepPointYZ_qpLast->Fill(y1,z1);
	    
	    h_timeOfLastStep->Fill(t1);
	  }	  
	}
      }
    }
  }  
  outFile->Write();
}




//----------------------------------------------------------------------------
// Analysis of QP Diffusion (Traping locations)
// (input = Validation_QPDiffusionPlusTrapping)
// This is literally the same analysis as the previous one, but for trapping
void ValidationAnalysis_QPDiffusionPlusTrapping() {
  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_QPDiffusionPlusTrapping.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge
  //interface
  TFile * outFile =
    new TFile("Validation_QPDiffusionPlusTrapping.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -2.6;
  double maxZ_mm = -1.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 10;
  double minTime_ns = 0;
  double maxTime_ns = 3e8;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_qpLast =
    new TH2F("h_postStepPointXY_qpLast",
	     "h_postStepPointXY, Last QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpLast =
    new TH2F("h_postStepPointXZ_qpLast",
	     "h_postStepPointXZ, Last QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpLast =
    new TH2F("h_postStepPointYZ_qpLast",
	     "h_postStepPointYZ, Last QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_qpLastYVsQPEnergy =
    new TH2F("h_qpLastYVsQPEnergy",
	     "h_qpLastZVsQPEnergy, Last QP Step Y vs. energy; QP Energy [meV]; Last QP Step Y",
	     nStepsE,minE_meV,maxE_meV,
	     1000,-(barLength_mm*1.1)/2.0,(barLength_mm*1.1)/2.0);
  
  TH1F * h_timeOfLastStep =
    new TH1F("h_timeOfLastStep","h_timeOfLastStep; Time [ns]; QPs;",
	     10000,minTime_ns,maxTime_ns);
  
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);

	  //Only fill this if the step is the last one
	  if (iS == theTrack.stepVect.size()-1) {
	    h_postStepPointXY_qpLast->Fill(x1,y1);
	    h_postStepPointXZ_qpLast->Fill(x1,z1);
	    h_postStepPointYZ_qpLast->Fill(y1,z1);
	    h_qpLastYVsQPEnergy->Fill(preStepKinEnergy_eV*1000,y1);	    
	    h_timeOfLastStep->Fill(t1);
	  }	  
	}
      }
    }
  }
  outFile->Write();
}

//----------------------------------------------------------------------------
// Analysis of QP Diffusion With Multiple Gaps
// (input = Validation_QPDiffusionMultiGap)
void ValidationAnalysis_QPDiffusionMultiGap() {
  //Read in the events from the textfile
  std::vector<Event> eventVect =
    ReadInG4CMPStepTextFile("Validation_QPDiffusionMultiGap.txt");
  std::cout << "Done filling events." << std::endl;

  //Primary goal of this analysis function: compute transmission at the Si-Ge
  //interface
  TFile * outFile = new TFile("Validation_QPDiffusionMultiGap.root","RECREATE");
  
  //Histogram definitions. 
  double barWidth_mm = 2.0;
  double barLength_mm = 10.0;
  double minZ_mm = -3.6;
  double maxZ_mm = -2.4;
  double nStepsE = 40;
  double minE_meV = 0;
  double maxE_meV = 11;
  double minTime_ns = 0;
  double maxTime_ns = 3e8;
  
  TH2F * h_postStepPointXY_qpAll =
    new TH2F("h_postStepPointXY_qpAll",
	     "h_postStepPointXY, All QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpAll =
    new TH2F("h_postStepPointXZ_qpAll",
	     "h_postStepPointXZ, All QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpAll =
    new TH2F("h_postStepPointYZ_qpAll",
	     "h_postStepPointYZ, All QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_qpLast =
    new TH2F("h_postStepPointXY_qpLast",
	     "h_postStepPointXY, Last QP Steps; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_qpLast =
    new TH2F("h_postStepPointXZ_qpLast",
	     "h_postStepPointXZ, Last QP Steps; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_qpLast =
    new TH2F("h_postStepPointYZ_qpLast",
	     "h_postStepPointYZ, Last QP Steps; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_phonons =
    new TH2F("h_postStepPointXY_phonons",
	     "h_postStepPointXY_phonons; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_phonons =
    new TH2F("h_postStepPointXZ_phonons",
	     "h_postStepPointXZ_phonons; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointYZ_phonons =
    new TH2F("h_postStepPointYZ_phonons",
	     "h_postStepPointYZ_phonons; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);

  TH2F * h_postStepPointXY_phononsLast =
    new TH2F("h_postStepPointXY_phononsLast",
	     "h_postStepPointXY_phononsLast; Step X [mm]; Step Y [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,-barLength_mm*1.1,barLength_mm*1.1);
  
  TH2F * h_postStepPointXZ_phononsLast =
    new TH2F("h_postStepPointXZ_phononsLast",
	     "h_postStepPointXZ_phononsLast; Step X [mm]; Step Z [mm]",
	     100,-barWidth_mm*1.1,barWidth_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
  
  TH2F * h_postStepPointYZ_phononsLast =
    new TH2F("h_postStepPointYZ_phononsLast",
	     "h_postStepPointYZ_phononsLast; Step Y [mm]; Step Z [mm]",
	     100,-barLength_mm*1.1,barLength_mm*1.1,
	     1000,minZ_mm,maxZ_mm);
    
  TH1F * h_timeOfLastStep_qp =
    new TH1F("h_timeOfLastStep_qp","h_timeOfLastStep_qp; Time [ns]; QPs;",
	     10000,minTime_ns,maxTime_ns);
  
  TH1F * h_energyOfLastStep_qp =
    new TH1F("h_energyOfLastStep_qp","h_energyOfLastStep_qp; Time [ns]; QPs;",
	     1000,minE_meV,maxE_meV);
  
  //Loop over events
  for (int iE = 0; iE < eventVect.size(); ++iE) {
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for (int iT = 0; iT < theEvent.trackVect.size(); ++iT) {
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over steps within tracks
      double lastRadStepTime_ns = 0;
      for (int iS = 0; iS < theTrack.stepVect.size(); ++iS) {
	Step theStep1 = theTrack.stepVect[iS];
       	
	//Compute quantities
	std::string processName = theStep1.processName;
	std::string particleName = theStep1.particleName;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double z_dir1 = dz_1 / fabs(dz_1);
	double preStepEnergy_eV = theStep1.preStepEnergy_eV;
	double preStepKinEnergy_eV = theStep1.preStepKinEnergy_eV;
	double postStepEnergy_eV = theStep1.postStepEnergy_eV;
	double postStepKinEnergy_eV = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	       	
	//Select only QPs
	if (particleName.find("QP") != std::string::npos) {

	  //Fill all-QP-step histograms. Note these are not physically
	  //meaningful because of the walk-on-spheres algorithm -- only steps
	  //with other physics processes (boundary hits,
	  //radiation, etc.) have physically meaningful positions.
	  h_postStepPointXY_qpAll->Fill(x1,y1);
	  h_postStepPointXZ_qpAll->Fill(x1,z1);
	  h_postStepPointYZ_qpAll->Fill(y1,z1);

	  //Only fill this if the step is the last one
	  if (iS == theTrack.stepVect.size()-1) {
	    h_postStepPointXY_qpLast->Fill(x1,y1);
	    h_postStepPointXZ_qpLast->Fill(x1,z1);
	    h_postStepPointYZ_qpLast->Fill(y1,z1);
	    h_timeOfLastStep_qp->Fill(t1);
	    h_energyOfLastStep_qp->Fill(preStepKinEnergy_eV*1000);
	  }	  
	}
	//Select only phonons
	if (particleName.find("phonon") != std::string::npos) {

	  //Fill all-phonon-step histograms.
	  h_postStepPointXY_phonons->Fill(x1,y1);
	  h_postStepPointXZ_phonons->Fill(x1,z1);
	  h_postStepPointYZ_phonons->Fill(y1,z1);

	  //Only fill this if the step is the last one
	  if (iS == theTrack.stepVect.size()-1) {
	    h_postStepPointXY_phononsLast->Fill(x1,y1);
	    h_postStepPointXZ_phononsLast->Fill(x1,z1);
	    h_postStepPointYZ_phononsLast->Fill(y1,z1);
	  }	  
	}
      }
    }
  }
  outFile->Write();
}

//-----------------------------------------------------------------------------
// Main function
void RunCompleteValidationAnalysis() {
  ValidationAnalysis_BoundaryTransmission();
  ValidationAnalysis_PolycrystalElasticScattering();
  ValidationAnalysis_SCPairbreaking();
  ValidationAnalysis_QPPhononRadiation();
  ValidationAnalysis_QPRecombination();
  ValidationAnalysis_QPDiffusionTOFP();
  ValidationAnalysis_QPDiffusionPlusRecombination();
  ValidationAnalysis_QPDiffusionPlusTrapping();
  ValidationAnalysis_QPDiffusionMultiGap();
}
