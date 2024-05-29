//---------------------------------------------------------
//
// RISQTutorialAnalysis.cc
// linehan3@fnal.gov
// 5/23/2024
//
// This is an analysis macro that is built to analyze
// the output of the G4CMP simulations that are run in
// the RISQ G4CMP Tutorial code.
//
//---------------------------------------------------------

//C++ includes
#include <cstdlib>
#include <fstream>
#include <iostream>

//ROOT includes
#include "TH2F.h"
#include "TH1F.h"

//---------------------------------------------------------------------------------------
// Define a set of structs for use interpreting the output from G4CMP
struct Hit
{
  int runID;
  int eventID;
  int trackID;
  std::string particleName;
  double startEnergy_eV;
  double startX_mm;
  double startY_mm;
  double startZ_mm;
  double startT_ns;
  double eDep_eV;
  double trackWeight;
  double endX_mm;
  double endY_mm;
  double endZ_mm;
  double endT_ns;
};

struct PrimaryInfo
{
  int runID;
  int eventID;
  int trackID;
  std::string particleName;
  double energy_eV;
  double X_mm;
  double Y_mm;
  double Z_mm;
  double T_ns;  
};

struct Event
{
  int runID;
  int eventID;
  std::vector<Hit> hitVect;
  PrimaryInfo thePrim;
};

//Forward declarations
std::vector<Event> ReadInG4CMPPrimaryAndHitFiles(std::string hitTextFilename, std::string primaryTextFilename);
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename);
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename);
int FindClosestQubitID(double hitX_mm, double hitY_mm);
  
void AnalyzeMuonEvent(std::string primariesFilename, std::string hitsFilename,double scaleFactorEHPairs)
{
  //Read in the text files
  std::vector<Event> eventList = ReadInG4CMPPrimaryAndHitFiles(hitsFilename,primariesFilename);
  std::cout << "Done reading in events." << std::endl;

  //Define an outfile
  TFile * fOut = new TFile("AnalysisOutput.root","RECREATE");
  
  //Define a number of histograms for the hits
  TH1F * h_eDep = new TH1F("h_eDep","Hit EDeps; log10(eDep[eV]); nEvents",200,-6,1);
  TH1F * h_qubitTotalHitEnergy_singleEvent = new TH1F("h_qubitTotalHitEnergy_singleEvent","In-Qubit Total Energy; Qubit ID; Energy [eV]",6,0,6);
  TH2F * h_hitXY = new TH2F("h_nHitsXY","XY Locations of Hits; X [mm]; Y [mm]; NHits/bin",200,-5,5,200,-5,5);
  TH2F * h_hitYZ = new TH2F("h_nHitsYZ","YZ Locations of Hits; Y [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_hitXZ = new TH2F("h_nHitsXZ","XZ Locations of Hits; X [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  
  //Loop over events
  for( int iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;
    
    //Get the event
    Event tE = eventList[iE];
    
    //Plot a number of hit-related things: 
    for ( int iH = 0; iH < tE.hitVect.size(); ++iH ){

      //Gather hit information
      double hitX_mm = tE.hitVect[iH].endX_mm;
      double hitY_mm = tE.hitVect[iH].endY_mm;
      double hitZ_mm = tE.hitVect[iH].endZ_mm;      
      double energy_eV = tE.hitVect[iH].eDep_eV;

      //From hit XY, figure out which qubit we're closest to, and add the
      //energy to that qubit id
      int closestQubitID = FindClosestQubitID(hitX_mm,hitY_mm);
      h_qubitTotalHitEnergy_singleEvent->Fill(closestQubitID,energy_eV);
      
      //Plot other hit information
      h_hitXY->Fill(hitX_mm,hitY_mm);
      h_hitYZ->Fill(hitY_mm,hitZ_mm);
      h_hitXZ->Fill(hitX_mm,hitZ_mm);
    }
  }

  //Last up, we need to remember that we did downsampling, so we need to scale our simulations back up to match
  h_qubitTotalHitEnergy_singleEvent->Scale(1.0/scaleFactorEHPairs);
  
  fOut->Write();
}

//---------------------------------------------------------------------------------------
// This is a bit of a kludge. What we really should do is pass the hit VOLUME out of
// G4CMP in the hit info. However, this involves some broader modifications to G4CMP
// beyond this example so we're kludging it together with this function.
int FindClosestQubitID(double hitX_mm, double hitY_mm)
{
  const int nQubits = 6;
  double qubitXY[nQubits][2] = {{-1.85,1.8},
				{0,1.8},
				{1.85,1.8},
				{-1.85,-1.8},
				{0,-1.8},
				{1.85,-1.8}};
  double lowestDeltaR2_mm = 1000000;
  int theQubitID = -1;
  for( int iQ = 0; iQ < nQubits; ++iQ ){
    double deltaR2_mm = TMath::Power((qubitXY[iQ][0] - hitX_mm),2) + TMath::Power((qubitXY[iQ][1] - hitY_mm),2);
    if( deltaR2_mm < lowestDeltaR2_mm ){
      lowestDeltaR2_mm = deltaR2_mm;
      theQubitID = iQ;
    }
  }
  return theQubitID;	 
}


//---------------------------------------------------------------------------------------
// Analysis Script: Phonon Collection Efficiency
void PCEStudy(std::string primariesFilename, std::string hitsFilename)
{
  //Read in the text files
  std::vector<Event> eventList = ReadInG4CMPPrimaryAndHitFiles(hitsFilename,primariesFilename);
  std::cout << "Done reading in events." << std::endl;


  //Define an outfile
  TFile * fOut = new TFile("AnalysisOutput.root","RECREATE");
  
  //Define a number of histograms for the hits
  TH1F * h_nHits = new TH1F("h_nHits","Number of Hits Per Event; nHits; nEvents",100,0,100);
  TH1F * h_eDep = new TH1F("h_eDep","Hit EDeps; log10(eDep[eV]); nEvents",200,-6,1);  
  TH2F * h_hitXY = new TH2F("h_nHitsXY","XY Locations of Hits; X [mm]; Y [mm]; NHits/bin",200,-5,5,200,-5,5);
  TH2F * h_hitYZ = new TH2F("h_nHitsYZ","YZ Locations of Hits; Y [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_hitXZ = new TH2F("h_nHitsXZ","XZ Locations of Hits; X [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  
  //Define a number of histograms for the primaries
  TH2F * h_primXY = new TH2F("h_primXY","XY Location of Primary; X [mm]; Y [mm]; Primaries/bin",200,-5,5,200,-5,5);
  TH2F * h_primYZ = new TH2F("h_primYZ","YZ Location of Primary; Y [mm]; Z [mm]; Primaries/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_primXZ = new TH2F("h_primXZ","XZ Location of Primary; X [mm]; Z [mm]; Primaries/bin",200,-5,5,200,4.6,5.02);

  //Define histograms useful for calculating PCE
  int nPCEBinsX = 50;
  int nPCEBinsY = 50;
  TH2F * h_totalHitEnergyAtPrimaryXY = new TH2F("h_totalHitEnergyAtPrimaryXY","Total Energy of Hits Generated by Primaries at this XY; X [mm]; Y [mm]; Energy/bin [eV]",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  TH2F * h_totalPrimaryEnergyAtPrimaryXY = new TH2F("h_totalPrimaryEnergyAtPrimaryXY","Total Energy of Primaries at this XY; X [mm]; Y [mm]; Energy/bin [eV]",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  TH2F * h_pceVsXY = new TH2F("h_pceVsXY","Phonon Collection Efficiency vs. Primary XY; X [mm]; Y [mm]; PCE",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  
  
  //Loop over events
  for( int iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;

    //Get the event
    Event tE = eventList[iE];
    
    //Add to the primary vector
    PrimaryInfo thePrim = tE.thePrim;
    h_primXY->Fill(thePrim.X_mm,thePrim.Y_mm);
    h_primXZ->Fill(thePrim.X_mm,thePrim.Z_mm);
    h_primYZ->Fill(thePrim.Y_mm,thePrim.Z_mm);
    h_totalPrimaryEnergyAtPrimaryXY->Fill(thePrim.X_mm,thePrim.Y_mm,thePrim.energy_eV);      
    
    //Plot a number of hit-related things: hit multiplicity, hit locations in XYZ, hits in XYZ weighted by energy, etc.
    h_nHits->Fill(tE.hitVect.size());
    for ( int iH = 0; iH < tE.hitVect.size(); ++iH ){

      //Gather hit information
      double hitX = tE.hitVect[iH].endX_mm;
      double hitY = tE.hitVect[iH].endY_mm;
      double hitZ = tE.hitVect[iH].endZ_mm;
      double startX = tE.hitVect[iH].startX_mm;
      double startY = tE.hitVect[iH].startY_mm;
      double startZ = tE.hitVect[iH].startZ_mm;
      double logEnergy_eV = TMath::Log10(tE.hitVect[iH].eDep_eV);


      
      //Plot hit information
      h_hitXY->Fill(hitX,hitY);
      h_hitYZ->Fill(hitY,hitZ);
      h_hitXZ->Fill(hitX,hitZ);
      h_eDep->Fill(logEnergy_eV);
      h_totalHitEnergyAtPrimaryXY->Fill(thePrim.X_mm,thePrim.Y_mm,tE.hitVect[iH].eDep_eV);
    }
  }

  //Post-processing division
  for( int iBX = 1; iBX <= h_totalHitEnergyAtPrimaryXY->GetNbinsX(); ++iBX ){
    for( int iBY = 1; iBY <= h_totalHitEnergyAtPrimaryXY->GetNbinsY(); ++iBY ){
      double num = h_totalHitEnergyAtPrimaryXY->GetBinContent(iBX,iBY);
      double denom = h_totalPrimaryEnergyAtPrimaryXY->GetBinContent(iBX,iBY);
      double pce = 0;
      if( denom != 0 ) pce = num/denom;
      h_pceVsXY->SetBinContent(iBX,iBY,pce);
    }
  }

  
  //Write stuff currently established
  fOut->Write();

}

  
  




  









//---------------------------------------------------------------------------------------
// Parsing function
std::vector<Event> ReadInG4CMPPrimaryAndHitFiles(std::string hitTextFilename, std::string primaryTextFilename)
{
  //Output: a vector of event objects nicely organizing our data.
  std::vector<Event> output;

  //First, let's open up our primary file and parse it. We'll get a map of int (event ID)
  //to primary info and a map of int (eventID) to a list of hits
  std::map<int,PrimaryInfo> primaryInfo = ParsePrimaryTextFileForPrimaries(primaryTextFilename);
  std::map<int,std::vector<Hit> > hitInfo = ParseHitTextFileForHits(hitTextFilename);

  //Now we do a final loop over event ID to merge these into actual events.
  for( std::map<int,PrimaryInfo>::iterator it = primaryInfo.begin(); it != primaryInfo.end(); ++it ){
    Event theEvent;
    theEvent.runID = it->second.runID;
    theEvent.eventID = it->first;
    theEvent.hitVect = hitInfo[it->first];
    theEvent.thePrim = it->second;
    output.push_back(theEvent);
  }
  return output;

}
  

//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename)
{
  std::map<int,std::vector<Hit> > output;
  std::vector<Hit> dummy;
  
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  int eventID = -1;
  int runID = -1;
  int counter = 0;
  while(1){
    if(!infile.good()) break;
    if(infile.is_open()){
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between spaces)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while(getline(check1,token,' ')){
	tokens.push_back(token);
      }
      if( tokens.size() == 0 ) break;
	  
      //If we're on the first line of the file
      if( tokens[0].find("Run") != std::string::npos ){
	continue;
      }

      //Check the runID and eventID, and if different than existing one,
      //push back a new event into the map
      if( std::atoi(tokens[0].c_str()) != runID || std::atoi(tokens[1].c_str()) != eventID ){
	output.emplace(std::atoi(tokens[1].c_str()),dummy);
	runID = std::atoi(tokens[0].c_str());
	eventID = std::atoi(tokens[1].c_str());
	counter++; 
	if( counter % 1000 == 0 ) std::cout << "Done reading " << counter << " events for hits." << std::endl;
      }

      //Log the hit information and push back into the most recently-created event in the vector
      Hit theHit;      
      theHit.runID = std::atoi(tokens[0].c_str());
      theHit.eventID = std::atoi(tokens[1].c_str());
      theHit.trackID = std::atoi(tokens[2].c_str());
      theHit.particleName = tokens[3];
      theHit.startEnergy_eV = std::atof(tokens[4].c_str());
      theHit.startX_mm = std::atof(tokens[5].c_str());
      theHit.startY_mm = std::atof(tokens[6].c_str());
      theHit.startZ_mm = std::atof(tokens[7].c_str());
      theHit.startT_ns = std::atof(tokens[8].c_str());
      theHit.eDep_eV = std::atof(tokens[9].c_str());
      theHit.trackWeight = std::atof(tokens[10].c_str());
      theHit.endX_mm = std::atof(tokens[11].c_str());
      theHit.endY_mm = std::atof(tokens[12].c_str());
      theHit.endZ_mm = std::atof(tokens[13].c_str());
      theHit.endT_ns = std::atof(tokens[14].c_str());
      output[eventID].push_back(theHit);
    }
  }
  return output;
}


//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename)
{
  std::map<int,PrimaryInfo> output;
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  while(1){
    if(!infile.good()) break;
    if(infile.is_open()){
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between spaces)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while(getline(check1,token,' ')){
	tokens.push_back(token);
      }
      if( tokens.size() == 0 ) break;
	  
      //If we're on the first line of the file
      if( tokens[0].find("Run") != std::string::npos ){
	continue;
      }

      //Here, it's simpler since there's one line per event (assuming only a single run)
      //So we can use the event id as an index.
      PrimaryInfo thePrim;
      thePrim.runID = std::atoi(tokens[0].c_str());
      thePrim.eventID = std::atoi(tokens[1].c_str());
      thePrim.particleName = tokens[2];
      thePrim.energy_eV = std::atof(tokens[3].c_str());
      thePrim.X_mm = std::atof(tokens[4].c_str());
      thePrim.Y_mm = std::atof(tokens[5].c_str());
      thePrim.Z_mm = std::atof(tokens[6].c_str());
      thePrim.T_ns = std::atof(tokens[7].c_str());
      output.emplace(thePrim.eventID,thePrim);
    }
  }
  return output;
}



