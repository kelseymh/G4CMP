////////////////////////////////////////////////////////////////////
//
// G4CMPStepConverter.cc
// linehan3@fnal.gov
//
// This script takes stepping output from running G4CMP
// and organizes it for easy analysis. Eventually we'll want
// to decouple the output parsing and the analysis, but
// that's mostly once we get going with high-stats sims.
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"


//---------------------------------------------------------------------------------------
// Define a set of structs for use interpreting the output from G4CMP
struct Step
{
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

struct Track
{
  int runID;
  int eventID;
  int trackID;
  std::vector<Step> stepVect;
};

struct Event
{
  int runID;
  int eventID;
  std::vector<Track> trackVect;
};



//Forward declarations
std::vector<Event> ReadInG4CMPStepTextFile(std::string filename);




//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that basic phonon
//transport is working as expected
void AnalyzeConvertedSteps_PhononBoundaryProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  
  //Useful definitions
  double diskRadius_mm = 38.1;
  
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_PhononBoundaryProcessStudy.root","RECREATE");
  
  //Plotting histos: general
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  
  //Plotting histos: internal boundaries
  TH1F * h_stepLengthAfterLongStep = new TH1F("h_stepLengthAfterLongStep","Step Length After Long Step (internal boundary hitters); Length [mm]; NSteps;",10000,0,40.);
  TH1F * h_stepDzMinusLastStepDzNorm = new TH1F("h_stepDzMinusLastStepDzNorm","Sign of Step Dz - Sign of Last Step Dz; Value; NSteps;",100,-3.0,3.0);
  
  //Plotting histos: non-internal boundaries: we should look and see that wall steps are reflected.
  //reflected. Looking for the step length after long steps that don't fall into those cuts is good.
  TH1F * h_stepLengthAfterLongStep_Walls = new TH1F("h_stepLengthAfterLongStep_Walls","Step Length After Long Step (Wall-hitters); Length [mm]; NSteps;",10000,0.,40.);
  

  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over pairs of steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size()-1; ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	Step theStep2 = theTrack.stepVect[iS+1];

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
	double stepLength2_mm = TMath::Power(dX_2*dX_2 + dY_2*dY_2 + dZ_2*dZ_2,0.5 );
	double dz_2 = z1_2-z0_2;

	double Delta_Dz_norm = dz_2/fabs(dz_2) - dz_1/fabs(dz_1);
	
	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);
	
	//First conditional, to get a pure set of internal-boundary-hitters: if the step is long and the step end is within 0.9 of the radius
	//of the disk and not on one of the ends, plot the length of the next step
	if( stepLength1_mm > 0.1 ){
	  if( r1 < (0.99 * diskRadius_mm) && z1 > -12.6 && z1 < 20.8 ){
	    h_stepLengthAfterLongStep->Fill(stepLength2_mm);
	    h_stepDzMinusLastStepDzNorm->Fill(Delta_Dz_norm);
	    std::cout << "dz_2: " << dz_2 << ", dz_1: " << dz_1 << std::endl;
	  }
	}

	//Second conditional, to get a pure set of external-boundary-hitters: as long as the step is not near an internal boundary, it's game
	if( stepLength1_mm > 0.1 ){
	  if( !(fabs(z1-12.7)<0.01) && !(fabs(z1-20.7)<0.01) ){
	    h_stepLengthAfterLongStep_Walls->Fill(stepLength2_mm);
	    if( stepLength2_mm > 0.1 ) {
	      std::cout << "Event: " << theEvent.eventID << ", Second conditional triggered with next-step length != 0: " << stepLength2_mm << ". X1,Y1,Z1 of first track: " << x1 << "," << y1 << "," << z1 << std::endl;
	    }
	  }
	}
      }
    }
  }

  outFile->Write();
 						       
}



//----------------------------------------------------------------------------------------
//This is an analysis that does validation of phonon polycrystalline scattering in a superconductor/normal metal
void AnalyzeConvertedSteps_PhononPolycrystallineElasticScatteringStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;

  //Output file
  TFile * outFile = new TFile("SteppingStudies_PhononPolycrystallineElasticScatteringStudy.root","RECREATE");
  
  //Plotting histos: general
  double minX_mm = -3.0e-2;
  double maxX_mm = 3.0e-2;
  double minY_mm = -3.0e-2;
  double maxY_mm = 3.0e-2;
  double minZ_mm = 12.6999;
  double maxZ_mm = 12.7249;
  double nBinsX = 1000;
  double nBinsY = 1000;
  double nBinsZ = 1000;
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",nBinsX,minX_mm,maxX_mm,nBinsY,minY_mm,maxY_mm);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",nBinsX,minX_mm,maxX_mm,nBinsZ,minZ_mm,maxZ_mm);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",nBinsY,minY_mm,maxY_mm,nBinsZ,minZ_mm,maxZ_mm);

  //Plotting histos: step length and polar angle (from dot product)
  double minL_nm = 0;
  double maxL_nm = 600;
  double nStepsL = 1000;
  TH1F * h_stepLength = new TH1F("h_stepLength","Step Length; Length [nm]; Steps;",nStepsL,minL_nm,maxL_nm);
  TH1F * h_cosPolarAngleOfScatter = new TH1F("h_cosPolarAngleOfScatter","Cos(Polar angle of scatter); Cos(v1#bullet v2); Scatters",1000,-1.0,1.0);
  TH1F * h_deltaE = new TH1F("h_deltaE","Step #Delta E [eV]; #Delta E [eV];Steps",1000,-0.005,0.005);
  TH1F * h_interStepDeltaE = new TH1F("h_interStepDeltaE","h_interStepDeltaE; Step-to-step (pre-step) #Delta E [eV];Steps",1000,-0.005,0.005);
  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      
      //Loop over pairs of steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size()-1; ++iS ){
	Step theStep = theTrack.stepVect[iS];
	Step theStep2 = theTrack.stepVect[iS+1];
	
	//Compute quantities
	double x0 = theStep.preStepX_mm;
	double y0 = theStep.preStepY_mm;
	double z0 = theStep.preStepZ_mm;
	double x1 = theStep.postStepX_mm;
	double y1 = theStep.postStepY_mm;
	double z1 = theStep.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double e0 = theStep.preStepEnergy_eV;
	double e1 = theStep.postStepEnergy_eV;
	double deltaE_eV = e1-e0;

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
	double stepLength2_mm = TMath::Power(dX_2*dX_2 + dY_2*dY_2 + dZ_2*dZ_2,0.5 );
	double dz_2 = z1_2-z0_2;
	double e0_2 = theStep2.preStepEnergy_eV;
	double interStepDeltaE_eV = e0_2 -e0;

	if(interStepDeltaE_eV != 0.0 ){
	  std::cout << "iS: " << iS << ", Inter-step deltaE = " << interStepDeltaE_eV << std::endl;
	}
	
	//Compute dot product
	double dotProductNormalized = (dX*dX_2 + dY*dY_2 + dZ*dZ_2)/stepLength2_mm/stepLength1_mm;
	
	
	
	
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);

	h_stepLength->Fill(stepLength1_mm*1.0e6);
	h_cosPolarAngleOfScatter->Fill(dotProductNormalized);
	h_deltaE->Fill(deltaE_eV);
	h_interStepDeltaE->Fill(interStepDeltaE_eV);
	
      }
    }
  }

  //Do some fits to extract the attenuation length
  TF1 * exponentialFit = new TF1("exponential","[0]*TMath::Exp(-x/[1])",minL_nm,maxL_nm);
  exponentialFit->SetParameter(0,1000);
  exponentialFit->SetParameter(1,30);
  h_stepLength->Fit(exponentialFit,"","",minL_nm,maxL_nm);
  outFile->Write();
}
  
//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that basic QP
//transport is working as expected
void AnalyzeConvertedSteps_BogoliubovQPBoundaryProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;

  //Useful definitions
  double diskRadius_mm = 38.1;
  double smallDiskRadius = 10.0;
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPBoundaryProcessStudy.root","RECREATE");

  //Plotting histos: general
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  
  //Plotting histos: internal boundaries
  TH1F * h_stepLengthAfterLongStep = new TH1F("h_stepLengthAfterLongStep","Step Length After Long Step (internal boundary hitters); Length [mm]; NSteps;",10000,0,40.);
  TH1F * h_stepDzMinusLastStepDzNorm = new TH1F("h_stepDzMinusLastStepDzNorm","Sign of Step Dz - Sign of Last Step Dz; Value; NSteps;",100,-3.0,3.0);
  TH1F * h_trackTotalNBounces = new TH1F("h_trackNBounces","Max nBounces in a track; nBounces; NTracks",30,0,30);
  TH1F * h_stepLengthAfterLongStep_halfCylinder = new TH1F("h_stepLengthAfterLongStep_halfCylinder","Step Length After Long Step (half-cylinder internal boundary hitters); Length [mm]; NSteps;",10000,0,80.);
  TH1F * h_stepDzMinusLastStepDzNorm_halfCylinder = new TH1F("h_stepDzMinusLastStepDzNorm_halfCylinder","Sign of Step Dz - Sign of Last Step Dz, half-Cylinder internal hitters; Value; NSteps;",100,-3.0,3.0);
    TH1F * h_stepLengthAfterLongStep_internalCylinder = new TH1F("h_stepLengthAfterLongStep_internalCylinder","Step Length After Long Step (internal-cylinder internal boundary hitters); Length [mm]; NSteps;",10000,0,80.);
  TH1F * h_stepDzMinusLastStepDzNorm_internalCylinder = new TH1F("h_stepDzMinusLastStepDzNorm_internalCylinder","Sign of Step Dz - Sign of Last Step Dz, internal-Cylinder internal hitters; Value; NSteps;",100,-3.0,3.0);


  
  
  //Plotting histos: non-internal boundaries: we should look and see that wall steps are reflected.
  //reflected. Looking for the step length after long steps that don't fall into those cuts is good.
  TH1F * h_stepLengthAfterLongStep_Walls = new TH1F("h_stepLengthAfterLongStep_Walls","Step Length After Long Step (Wall-hitters); Length [mm]; NSteps;",10000,0.,40.);

  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over pairs of steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size()-1; ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	Step theStep2 = theTrack.stepVect[iS+1];

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
	double stepLength2_mm = TMath::Power(dX_2*dX_2 + dY_2*dY_2 + dZ_2*dZ_2,0.5 );
	double dz_2 = z1_2-z0_2;
	double nBounces = theStep2.nReflections;	
	
	double Delta_Dz_norm = dz_2/fabs(dz_2) - dz_1/fabs(dz_1);
	
	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);

	if( nBounces > maxNBounces ){maxNBounces = nBounces; }	
	
	//First conditional, to get a pure set of internal-boundary-hitters: if the step is long and the step end is within 0.9 of the radius
	//of the disk and not on one of the ends, plot the length of the next step
	if( stepLength1_mm > 0.1 ){
	  if( r1 < (0.99 * diskRadius_mm) && z1 > -12.6 && z1 < 20.8 ){
	    h_stepLengthAfterLongStep->Fill(stepLength2_mm);
	    h_stepDzMinusLastStepDzNorm->Fill(Delta_Dz_norm);
	    std::cout << "dz_2: " << dz_2 << ", dz_1: " << dz_1 << ", dz_2/fabs(dz_2): " << dz_2/fabs(dz_2) << std::endl;
	  }
	}

	//Second conditional, to get a pure set of internal-boundary hitters hitting the boundary between two half-cylinders. If the step is long
	//and the end is within 0.9 of the radius of the disk, not near a z boundary, and within a tolerance of y=0, then plot the length of the
	//next step.
	if( stepLength1_mm > 0.1 ){
	  if( r1 < (0.99 * diskRadius_mm) && z1 > 12.8 && z1 < 20.6 && fabs(y1) < 0.01 ){
	    h_stepLengthAfterLongStep_halfCylinder->Fill(stepLength2_mm);
	    h_stepDzMinusLastStepDzNorm_halfCylinder->Fill(Delta_Dz_norm);
	  }
	}

	//Third conditional: pure set of internal-boundary hitters hitting the boundary between a cylinder and an internal cylinder
	if( stepLength1_mm > 0.1 ){
	  if( r1 < (smallDiskRadius * 1.01) && z1 > 12.8 && z1 < 20.6 ){
	    h_stepLengthAfterLongStep_internalCylinder->Fill(stepLength2_mm);
	    h_stepDzMinusLastStepDzNorm_internalCylinder->Fill(Delta_Dz_norm);
	  }
	}
	
	//Second conditional, to get a pure set of external-boundary-hitters: as long as the step is not near an internal boundary, it's game
	if( stepLength1_mm > 0.1 ){
	  if( !(fabs(z1-12.7)<0.01) && !(fabs(z1-20.7)<0.01) ){
	    h_stepLengthAfterLongStep_Walls->Fill(stepLength2_mm);
	    if( stepLength2_mm > 0.1 ) {
	      std::cout << "Event: " << theEvent.eventID << ", Second conditional triggered with next-step length != 0: " << stepLength2_mm << ". X1,Y1,Z1 of first track: " << x1 << "," << y1 << "," << z1 << std::endl;
	    }
	  }
	}
      }
      h_trackTotalNBounces->Fill(maxNBounces);
    }
  }

  outFile->Write();
 						       
}



//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that basic QP
//transport is working as expected
void AnalyzeConvertedSteps_BogoliubovQPLocalTrappingProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In Local Trapping Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;
  
  //Useful definitions
  double diskRadius_mm = 38.1;
  double minL_mm = 0;
  double maxL_mm = 10;
  double nStepsL = 100;
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPLocalTrappingProcessStudy.root","RECREATE");

  //Plotting histos: general
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_stepLength = new TH1F("h_stepLength","Step Length; Length [mm]; Steps;",nStepsL,minL_mm,maxL_mm);
  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
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

	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);
	h_stepLength->Fill(stepLength1_mm);
      }
    }
  }



  //Do some fits to extract the attenuation length
  TF1 * exponentialFit = new TF1("exponential","[0]*TMath::Exp(-x/[1])",minL_mm,maxL_mm);
  exponentialFit->SetParameter(0,10);
  exponentialFit->SetParameter(1,1);
  h_stepLength->Fit(exponentialFit,"","",minL_mm,maxL_mm);
  outFile->Write();
 						       
}

//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that basic QP
//transport is working as expected
void AnalyzeConvertedSteps_BogoliubovQPLocalTrappingWithDiffusionProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In Local Trapping Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;
  
  //Useful definitions
  double diskRadius_mm = 38.1;
  double minL_mm = 0;
  double maxL_mm = 10;
  double nStepsL = 100;
  double minT_ns = 0;
  double maxT_ns = 6e6;
  double nStepsT = 1000;

  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPLocalTrappingWithDiffusionProcessStudy.root","RECREATE");

  //Plotting histos: general
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_stepLength = new TH1F("h_stepLength","Step Length; Length [mm]; Steps;",nStepsL,minL_mm,maxL_mm);
  TH1F * h_stepDuration = new TH1F("h_stepDuration","Step Duration; Duration [ns]; Steps;",nStepsT,minT_ns,maxT_ns);
  TH1F * h_deathTime = new TH1F("h_deathTime","Death Time; Time [ns]; Steps;",nStepsT,minT_ns,maxT_ns);
  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];

	//Compute quantities
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double dz_1 = z1-z0;
	double dur_ns = t1-t0;

	//General plots
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);
	h_stepLength->Fill(stepLength1_mm);
	h_stepDuration->Fill(dur_ns);

	//Calculate death time
	if(iS == theTrack.stepVect.size()-1 ){
	  std::cout << "deathTime: " << t1 << std::endl;
	  h_deathTime->Fill(t1);
	}
      }
    }
  }



  //Do some fits to extract the attenuation length
  TF1 * exponentialFit = new TF1("exponential","[0]*TMath::Exp(-x/[1])",minT_ns,maxT_ns);
  exponentialFit->SetParameter(0,10);
  exponentialFit->SetParameter(1,1e6);
  h_deathTime->Fit(exponentialFit,"","",minT_ns,maxT_ns);
  outFile->Write();
 						       
}



//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that basic phonon
//transport is working as expected
void AnalyzeConvertedSteps_PhononSCPairbreakingProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In Phonon SC Pairbreaking Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;

  double gap_meV = 0.176; //Useful 
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_PhononSCPairbreakingProcessStudy.root","RECREATE");  

  //Defining variables
  int nStepsE = 20;
  double minE_meV = 0.0;
  double maxE_meV = 11;
  int nStepsDeltaT = 200;
  double minDeltaT_ns = 0;
  double maxDeltaT_ns = 1.0;
  double tau0_ph_ns = 0.242;
  int nStepsDist = 1000;
  double minDist_mm = 0;
  double maxDist_mm = 10;
  
  //For the step length, we'll want to rescale to get tau0/deltaT (where deltaT is step length divided by the phonon velocity)
  int nStepsTau0DivDeltaT = 200;
  double minTau0DivDeltaT = 0.0;
  double maxTau0DivDeltaT = 20.0;  
  
  //Defining histograms
  TH2F * h_postStepPointXY = new TH2F("h_postStepPointXY","PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_postStepPointXZ = new TH2F("h_postStepPointXZ","PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_postStepPointYZ = new TH2F("h_postStepPointYZ","PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_Tau0DivDeltaTVsPreStepEnergy = new TH2F("h_Tau0DivDeltaTVsPreStepEnergy","Tau0DivDeltaT vs Pre-Step Energy; Pre-Step Energy [meV]; Tau0DivDeltaT [mm]",nStepsE,minE_meV,maxE_meV,nStepsTau0DivDeltaT,minTau0DivDeltaT,maxTau0DivDeltaT);
  TH2F * h_DeltaTVsPreStepEnergy = new TH2F("h_DeltaTVsPreStepEnergy","DeltaT vs Pre-Step Energy; Pre-Step Energy [meV]; DeltaT [mm]",nStepsE,minE_meV,maxE_meV,nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH1F * h_deltaT = new TH1F("h_deltaT","Step DeltaT; DeltaT [ns]; Steps",nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH1F * h_preStepEnergy = new TH1F("h_preStepEnergy","Pre-step energy; Energy [meV]; Steps",nStepsE,minE_meV,maxE_meV);
  TH1F * h_postStepEnergy = new TH1F("h_postStepEnergy","Post-step energy; Energy [meV]; Steps",nStepsE,minE_meV,maxE_meV);

  //Histograms for QPs
  TH2F * h_qpEnergyVsPreStepEnergy = new TH2F("h_qpEnergyVsPreStepEnergy","QP Energy vs. Pre-Step (Phonon) Energy; Phonon energy [meV]; QP Energies [meV]",nStepsE*10,minE_meV,maxE_meV,nStepsE*10,minE_meV,maxE_meV);
  TH1F * h_qpDistToProgenitorPhononEnd = new TH1F("h_qpDistToProgenitorPhononEnd","Distance between QP start and progenitor phonon end; Distance [mm]",nStepsDist,minDist_mm,maxDist_mm);
	  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    //Do a preliminary loop over tracks to find the phononTF and extract its energy
    double phononTFEnergy_meV = 0;
    double pLX = 0;
    double pLY = 0;
    double pLZ = 0;
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      bool isPhononTF = false;
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	std::string particleName = theStep1.particleName;
	if( particleName == "phononTF" ){
	  phononTFEnergy_meV = theStep1.preStepKinEnergy_eV * 1000;
	  isPhononTF = true;
	}
	//If we're on the last step, take the final track point
	if( isPhononTF && iS == theTrack.stepVect.size()-1 ){
	  pLX = theStep1.postStepX_mm;
	  pLY = theStep1.postStepY_mm;
	  pLZ = theStep1.postStepZ_mm;
	}
      }
    }
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double e0 = theStep1.preStepKinEnergy_eV;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double e1 = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double stepDeltaT_ns = t1-t0;
	double dz_1 = z1-z0;
	std::string particleName = theStep1.particleName;

	//Compute the scaled rate
	double tau0DivDeltaT = tau0_ph_ns / stepDeltaT_ns;

	if( particleName == "BogoliubovQP" ){
	  h_qpEnergyVsPreStepEnergy->Fill(phononTFEnergy_meV,e0*1000);

	  //Compute the distance between the starting point of this QP and the end point of the progenitor phonon
	  double deltaX_mm = x0 - pLX;
	  double deltaY_mm = y0 - pLY;
	  double deltaZ_mm = z0 - pLZ;
	  double dist_mm = TMath::Power(deltaX_mm*deltaX_mm + deltaY_mm * deltaY_mm + deltaZ_mm * deltaZ_mm,0.5);
	  h_qpDistToProgenitorPhononEnd->Fill(dist_mm);  
	}

	//Filter out tracks that are not phonons (there should be two QP tracks per event, but we will kill these for this study)
	if( particleName != "phononTF") continue;

	//First up, basic plots: starting/stopping points
	h_postStepPointXY->Fill(x1,y1);
	h_postStepPointXZ->Fill(x1,z1);
	h_postStepPointYZ->Fill(y1,z1);	
	h_deltaT->Fill(stepDeltaT_ns);
	
	//Now plot the step length vs kinetic energy of the particle in the pre-step point
	h_preStepEnergy->Fill(e0*1000); //To get into meV
	h_postStepEnergy->Fill(e1*1000); //To get into meV
	h_Tau0DivDeltaTVsPreStepEnergy->Fill(e0*1000,tau0DivDeltaT);
	h_DeltaTVsPreStepEnergy->Fill(e0*1000,stepDeltaT_ns);
      }
    }
  }

  //Set up fitter to do lots of iterations
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000);
  
  //Do fits to vertical slices of the deltaTVsPreStepEnergy plot, and get the decay time
  double fitMin_ns = 0.01;
  //double fitMax_ns = 0.4; //for low energy
  double fitMax_ns = 0.2; //For high energy
  TF1 * expFit = new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",fitMin_ns,fitMax_ns);
  TGraphErrors * g_fitExponents = new TGraphErrors();
  int offset = 0; //For plotting
  for( int iBX = 1; iBX < h_DeltaTVsPreStepEnergy->GetNbinsX(); ++iBX){
    std::cout << " pre-step energy bin center: " << h_DeltaTVsPreStepEnergy->GetXaxis()->GetBinCenter(iBX) << std::endl;
    if( h_DeltaTVsPreStepEnergy->GetXaxis()->GetBinCenter(iBX) < 2*gap_meV ){
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }
    
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice = (TH1F*)h_DeltaTVsPreStepEnergy->ProjectionY(name,iBX,iBX,"");
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");
    //expFit->SetParameter(0,2000./(theSlice->GetMean()*60)); //Hardcoded start point for low energy 
    //expFit->SetParameter(1,0.35); //Hardcoded start point for low-energy
    expFit->SetParameter(0,2000); //Hardcoded start point (based on number of entries)
    expFit->SetParameter(1,0.2); //Hardcoded start point for high-energy
    theSlice->Fit(expFit,"","",fitMin_ns,fitMax_ns);
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);

    std::cout << "Fit exponent: " << fitExponent_ns << std::endl;
    
    //Convert the fit exponent to tau0/tauB
    double tau0DivTauB = tau0_ph_ns / fitExponent_ns;
    double tau0DivTauB_error = tau0DivTauB * (fitExponentError_ns / fitExponent_ns);
    std::cout << "Tau0DivTauB = " << tau0DivTauB << ", and plotting." << std::endl;
    
    g_fitExponents->SetPoint(iBX-offset-1,h_DeltaTVsPreStepEnergy->GetXaxis()->GetBinCenter(iBX),tau0DivTauB);
    g_fitExponents->SetPointError(iBX-offset-1,0,tau0DivTauB_error);
  }
  g_fitExponents->SetName("g_fitExponents");
  g_fitExponents->SetTitle("Fit Exponents");
  g_fitExponents->GetXaxis()->SetTitle("Pre-Step Energy [meV]");
  g_fitExponents->GetYaxis()->SetTitle("Tau0 / TauB");

  TCanvas * c0 = new TCanvas();
  g_fitExponents->Draw("ALP");
  
  g_fitExponents->Write();
  outFile->Write();
}



//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that phonon radiation by
//quasiparticles is working as expected
void AnalyzeConvertedSteps_BogoliubovQPRadiatesPhononProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In BogoliubovQP Radiates Phonons Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;

  double gap_meV = 0.176; //Useful, gap0
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPRadiatesPhononProcessStudy.root","RECREATE");  

  //Defining variables
  int nStepsE = 20;
  double minE_meV = 0.0;
  double maxE_meV = 11;
  int nStepsR = 200;
  double minR_mm = 0.0;
  double maxR_mm = 40.0;
  int nStepsZ = 200;
  double minZ_mm = -12.8;
  double maxZ_mm = 25.0;
  int nStepsDeltaT = 5000;
  double minDeltaT_ns = 0;
  double maxDeltaT_ns = 500;
  int nStepsDist = 1000;
  double minDist_mm = 0;
  double maxDist_mm = 10;
  int nTheta = 200;
  double minTheta_deg = 0.;
  double maxTheta_deg = 180;
  int nPhi = 200;
  double minPhi_deg = -180.;
  double maxPhi_deg = 180.;
  double tau0_qp_ns = 438.;
  
  //Defining generic QP histograms (all steps)
  TH2F * h_QP_postStepPointXY = new TH2F("h_QP_postStepPointXY","QP PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QP_postStepPointXZ = new TH2F("h_QP_postStepPointXZ","QP PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QP_postStepPointYZ = new TH2F("h_QP_postStepPointYZ","QP PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_QP_deltaT = new TH1F("h_QP_deltaT","QP Step DeltaT; DeltaT [ns]; Steps",nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH1F * h_QP_preStepEnergy_step0 = new TH1F("h_QP_preStepEnergy_step0","QP Pre-step energy, step 0; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_postStepEnergy_step0 = new TH1F("h_QP_postStepEnergy_step0","QP Post-step energy, step 0; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_preStepEnergy = new TH1F("h_QP_preStepEnergy","QP Pre-step energy; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_postStepEnergy = new TH1F("h_QP_postStepEnergy","QP Post-step energy; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_finalStepPostStepEnergy = new TH1F("h_QP_finalStepPostStepEnergy","QP Final Step Post-step energy; Energy [meV]; Counts",nStepsE*400,minE_meV,maxE_meV);
  TH2F * h_QP_finalStepPostStepEnergyVsR = new TH2F("h_QP_finalStepPostStepEnergyVsR","QP Final Step Post-step energy vs. R; Energy [meV]; R [mm]; Counts",nStepsE*400,minE_meV,maxE_meV,nStepsR,0,maxR_mm);
  TH2F * h_QP_finalStepPostStepEnergyVsZ = new TH2F("h_QP_finalStepPostStepEnergyVsZ","QP Final Step Post-step energy vs. Z; Energy [meV]; Z [mm]; Counts",nStepsE*400,minE_meV,maxE_meV,nStepsZ,minZ_mm,maxZ_mm);
  
  
  
  //Defining QP histograms in steps where the step is limited by a BogoliubovQPRadiatesPhonon process
  TH2F * h_QPrad_postStepPointXY = new TH2F("h_QPrad_postStepPointXY","QP PostStepPoint XY, Phonon Radiation Steps; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QPrad_postStepPointXZ = new TH2F("h_QPrad_postStepPointXZ","QP PostStepPoint XZ, Phonon Radiation Steps; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QPrad_postStepPointYZ = new TH2F("h_QPrad_postStepPointYZ","QP PostStepPoint YZ, Phonon Radiation Steps; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QPrad_DeltaTVsPreStepEnergy = new TH2F("h_QPrad_DeltaTVsPreStepEnergy","QP DeltaT vs Pre-Step Energy, Phonon Radiation Steps; Pre-Step Energy [meV]; DeltaT [ns]",nStepsE*4,minE_meV,maxE_meV,nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH2F * h_QPrad_DeltaTVsPreStepEnergy_LowE = new TH2F("h_QPrad_DeltaTVsPreStepEnergy_LowE","QP DeltaT vs Pre-Step Energy, Phonon Radiation Steps; Pre-Step Energy [meV]; DeltaT [ns]",nStepsE*2,minE_meV,maxE_meV/5.,nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH1F * h_QPrad_deltaT = new TH1F("h_QPrad_deltaT","QP Step DeltaT, Phonon Radiation Steps; DeltaT [ns]; Steps",nStepsDeltaT,minDeltaT_ns,maxDeltaT_ns);
  TH1F * h_QPrad_preStepEnergy = new TH1F("h_QPrad_preStepEnergy","QP Pre-step energy, Phonon Radiation Steps; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QPrad_postStepEnergy = new TH1F("h_QPrad_postStepEnergy","QP Post-step energy, Phonon Radiation Steps; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  
  //Histograms for emergent phonons
  TH1F * h_phononPolarization = new TH1F("h_phononPolarization","Emitted Phonon Polarization; Value; Counts",3,0,3);
  TH2F * h_phononThetaPhi = new TH2F("h_phononThetaPhi","Emitted Phonon velocity theta and phi, Relative to Z axis; Theta [degrees]; Phi [degrees]",nTheta,minTheta_deg,maxTheta_deg,nPhi,minPhi_deg,maxPhi_deg);

  //Histograms for joing information about phonons and QPs
  TH2F * h_phononEnergyVsQPPreStepEnergy = new TH2F("h_phononEnergyVsQPPreStepEnergy","Phonon Energy vs. QP Pre-Step Energy; QP Pre-step energy [meV]; Phonon Energies [meV]",nStepsE*10,minE_meV,maxE_meV,nStepsE*10,minE_meV,maxE_meV);
  TH2F * h_corrPhononEnergyVsQPPreStepEnergy = new TH2F("h_corrPhononEnergyVsQPPreStepEnergy","(QPE-PE) vs. QP Pre-Step Energy; QP Pre-step energy [meV]; QP Energy - Phonon Energy [meV]",nStepsE*10,minE_meV,maxE_meV,nStepsE*10,minE_meV,maxE_meV);
  TH1F * h_changeInOverallEnergy = new TH1F("h_changeInOverallEnergy","Phonon+QP post-step energy - QP pre-step energy",10001,-maxE_meV/10.,maxE_meV/10.);
  
  //Okay, time to start looping over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    //Do a preliminary loop over tracks to find the initial BogoliubovQP and extract its energy
    double qpPreStepEnergy_meV = 0; //For just the first step
    double qpPostStepEnergy_meV = 0; //For just the first step
    double qpLX = 0;
    double qpLY = 0;
    double qpLZ = 0;
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      bool isBogoliubovQP = false;
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	std::string particleName = theStep1.particleName;
	if( particleName == "BogoliubovQP" && iS == 0 ){
	  qpPreStepEnergy_meV = theStep1.preStepKinEnergy_eV * 1000;
	  qpPostStepEnergy_meV = theStep1.postStepKinEnergy_eV * 1000;
	  isBogoliubovQP = true;
	}
	//If we're on the last step, take the final track point
	if( isBogoliubovQP && iS == theTrack.stepVect.size()-1 ){
	  qpLX = theStep1.postStepX_mm;
	  qpLY = theStep1.postStepY_mm;
	  qpLZ = theStep1.postStepZ_mm;
	}
      }
    }

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	int trackID = theStep1.trackID;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double e0 = theStep1.preStepKinEnergy_eV;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double e1 = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double dR = TMath::Power(x1*x1 + y1*y1 + z1*z1, 0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double stepDeltaT_ns = t1-t0;
	double dz_1 = z1-z0;
	std::string particleName = theStep1.particleName;
	std::string processName = theStep1.processName;

	//If we're looking at a QP, plot generic info
	if( particleName == "BogoliubovQP" ){
	  h_QP_postStepPointXY->Fill(x1,y1);
	  h_QP_postStepPointXZ->Fill(x1,z1);
	  h_QP_postStepPointYZ->Fill(y1,z1);
	  h_QP_deltaT->Fill(t1-t0);
	  h_QP_preStepEnergy->Fill(e0*1000);
	  h_QP_postStepEnergy->Fill(e1*1000);
	  if( iS == 0 ){
	    h_QP_preStepEnergy_step0->Fill(e0*1000.);
	    h_QP_postStepEnergy_step0->Fill(e1*1000.);
	  }
	  if( iS == theTrack.stepVect.size()-1 ){
	    h_QP_finalStepPostStepEnergy->Fill(e1*1000);
	    h_QP_finalStepPostStepEnergyVsR->Fill(e1*1000,r1);
	    h_QP_finalStepPostStepEnergyVsZ->Fill(e1*1000,z1);
	  }

	  //If we're looking at a QP *and* the step is also a qpradiatesphonon step is bogoliubovQPRadiatesPhonon, then a different
	  //set of histograms gets filled
	  if( processName == "bogoliubovQPRadiatesPhonon" ){
	    h_QPrad_postStepPointXY->Fill(x1,y1);
	    h_QPrad_postStepPointXZ->Fill(x1,z1);
	    h_QPrad_postStepPointYZ->Fill(y1,z1);
	    h_QPrad_DeltaTVsPreStepEnergy->Fill(e0*1000,t1-t0);
	    h_QPrad_DeltaTVsPreStepEnergy_LowE->Fill(e0*1000,t1-t0);
	    h_QPrad_deltaT->Fill(t1-t0);
	    h_QPrad_preStepEnergy->Fill(e0*1000);
	    h_QPrad_postStepEnergy->Fill(e1*1000);
	  }
	  
	}

	//If we're looking at phonons...
	if( particleName == "phononTF" || particleName == "phononTS" || particleName == "phononL" ){

	  //Associate the polarization to a number
	  double polariz = -1;
	  if( particleName == "phononL" ) polariz = 0.5;
	  if( particleName == "phononTF" ) polariz = 1.5;
	  if( particleName == "phononTS" ) polariz = 2.5;
	  h_phononPolarization->Fill(polariz);	  
	  h_phononThetaPhi->Fill(TMath::ACos(dZ/dR)*180./TMath::Pi(),TMath::ATan2(dY,dX)*180./TMath::Pi());
	  
	  //And, if we're looking at specifically the trackID = 2 particle (a phonon), then we can use the existing energy
	  //of the initial QP to get a distribution of which energies come from which other energies. Can't do this easily for
	  //later tracks because bookkeeping will be a nightmare
	  if( trackID == 2 ){
	    h_phononEnergyVsQPPreStepEnergy->Fill(qpPreStepEnergy_meV,e0*1000);
	    h_corrPhononEnergyVsQPPreStepEnergy->Fill(qpPreStepEnergy_meV,qpPreStepEnergy_meV-e0*1000);
	    h_changeInOverallEnergy->Fill(qpPreStepEnergy_meV-(qpPostStepEnergy_meV+e0*1000)); //Should always be zero here	    
	  }
	}
      }
    }
  }

  
  //Let's do some fits
  
  //More will go here in a bit
  //Set up fitter to do lots of iterations
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000);
  
  //Do fits to vertical slices of the deltaTVsPreStepEnergy plot, and get the decay time
  double fitMin_ns = 0.1;
  double fitMax_ns = 500; //For high energy
  TF1 * expFit = new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",fitMin_ns,fitMax_ns);
  TGraphErrors * g_fitExponents = new TGraphErrors();
  int offset = 0; //For plotting
  for( int iBX = 1; iBX < h_QPrad_DeltaTVsPreStepEnergy_LowE->GetNbinsX(); ++iBX){
    std::cout << " pre-step energy bin center: " << h_QPrad_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX) << std::endl;
    if( h_QPrad_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX) < 1*gap_meV ){
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }
    
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice = (TH1F*)h_QPrad_DeltaTVsPreStepEnergy_LowE->ProjectionY(name,iBX,iBX,"");
    if( h_QPrad_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX) < 4*gap_meV ){ theSlice->Rebin(5); }
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");
    //expFit->SetParameter(0,2000./(theSlice->GetMean()*60)); //Hardcoded start point for low energy 
    //expFit->SetParameter(1,0.35); //Hardcoded start point for low-energy
    expFit->SetParameter(0,theSlice->GetBinContent(1)); //Hardcoded start point (based on number of entries)
    expFit->SetParameter(1,theSlice->GetMean()); //Hardcoded start point for high-energy
    theSlice->Fit(expFit,"","",fitMin_ns,fitMax_ns);
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);

    std::cout << "Fit exponent: " << fitExponent_ns << std::endl;
    
    //Convert the fit exponent to taus/tau0
    double tauSDivTau0 = fitExponent_ns / tau0_qp_ns;
    double tauSDivTau0_error = tauSDivTau0 * (fitExponentError_ns / fitExponent_ns);
    std::cout << "TauSDivTau0 = " << tauSDivTau0 << ", and plotting." << std::endl;
    
    g_fitExponents->SetPoint(iBX-offset-1,h_QPrad_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX),tauSDivTau0);
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




//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that phonon radiation by
//quasiparticles is working as expected
void AnalyzeConvertedSteps_BogoliubovQPRecombinationProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In BogoliubovQP Recombination Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;

  double gap_meV = 0.176; //Useful, gap0
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPRecombinationProcessStudy.root","RECREATE");  

  //Defining variables
  int nStepsE = 20;
  double minE_meV = 0.0;
  double maxE_meV = 11;
  int nStepsR = 200;
  double minR_mm = 0.0;
  double maxR_mm = 40.0;
  int nStepsZ = 200;
  double minZ_mm = -12.8;
  double maxZ_mm = 25.0;
  int nStepsDeltaTLong = 5000;
  double minDeltaTLong_ns = 0;
  double maxDeltaTLong_ns = 10000; //500E4

  int nStepsDeltaT = 5000;
  double minDeltaT_ns = 0;
  double maxDeltaT_ns = 500;
  
  int nStepsDist = 1000;
  double minDist_mm = 0;
  double maxDist_mm = 10;
  int nTheta = 200;
  double minTheta_deg = 0.;
  double maxTheta_deg = 180;
  int nPhi = 200;
  double minPhi_deg = -180.;
  double maxPhi_deg = 180.;
  double tau0_qp_ns = 438.;
  
  //Defining generic QP histograms (all steps)
  TH2F * h_QP_postStepPointXY = new TH2F("h_QP_postStepPointXY","QP PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QP_postStepPointXZ = new TH2F("h_QP_postStepPointXZ","QP PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QP_postStepPointYZ = new TH2F("h_QP_postStepPointYZ","QP PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_QP_deltaT = new TH1F("h_QP_deltaT","QP Step DeltaT; DeltaT [ns]; Steps",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH1F * h_QP_preStepEnergy = new TH1F("h_QP_preStepEnergy","QP Pre-step energy; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_postStepEnergy = new TH1F("h_QP_postStepEnergy","QP Post-step energy; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  
  //Defining QP histograms in steps where the step is limited by a BogoliubovQPRadiatesPhonon process
  TH2F * h_QPrec_postStepPointXY = new TH2F("h_QPrec_postStepPointXY","QP PostStepPoint XY, Recombination Steps; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QPrec_postStepPointXZ = new TH2F("h_QPrec_postStepPointXZ","QP PostStepPoint XZ, Recombination Steps; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QPrec_postStepPointYZ = new TH2F("h_QPrec_postStepPointYZ","QP PostStepPoint YZ, Recombination Steps; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QPrec_DeltaTVsPreStepEnergy = new TH2F("h_QPrec_DeltaTVsPreStepEnergy","QP DeltaT vs Pre-Step Energy, Recombination Steps; Pre-Step Energy [meV]; DeltaT [ns]",nStepsE*4,minE_meV,maxE_meV,nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH2F * h_QPrec_DeltaTVsPreStepEnergy_LowE = new TH2F("h_QPrec_DeltaTVsPreStepEnergy_LowE","QP DeltaT vs Pre-Step Energy, Recombination Steps; Pre-Step Energy [meV]; DeltaT [ns]",nStepsE*2,minE_meV,maxE_meV/5.,nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH1F * h_QPrec_deltaT = new TH1F("h_QPrec_deltaT","QP Step DeltaT, Recombination Steps; DeltaT [ns]; Steps",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH1F * h_QPrec_preStepEnergy = new TH1F("h_QPrec_preStepEnergy","QP Pre-step energy, Recombination Steps; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QPrec_postStepEnergy = new TH1F("h_QPrec_postStepEnergy","QP Post-step energy, Recombination Steps; Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  
  //Histograms for emergent phonons
  TH1F * h_phononPolarization = new TH1F("h_phononPolarization","Emitted Phonon Polarization; Value; Counts",3,0,3);
  TH2F * h_phononThetaPhi = new TH2F("h_phononThetaPhi","Emitted Phonon velocity theta and phi, Relative to Z axis; Theta [degrees]; Phi [degrees]",nTheta,minTheta_deg,maxTheta_deg,nPhi,minPhi_deg,maxPhi_deg);

  //Histograms for joing information about phonons and QPs
  TH2F * h_phononEnergyVsQPPreStepEnergy = new TH2F("h_phononEnergyVsQPPreStepEnergy","Phonon Energy vs. QP Pre-Step Energy; QP Pre-step energy [meV]; Phonon Energies [meV]",nStepsE*10,minE_meV,maxE_meV,nStepsE*10,minE_meV,maxE_meV);
  
  //Okay, time to start looping over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    //Do a preliminary loop over tracks to find the initial BogoliubovQP and extract its energy
    double qpPreStepEnergy_meV = 0; //For just the first step
    double qpPostStepEnergy_meV = 0; //For just the first step
    double qpLX = 0;
    double qpLY = 0;
    double qpLZ = 0;
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      bool isBogoliubovQP = false;
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	std::string particleName = theStep1.particleName;
	if( particleName == "BogoliubovQP" && iS == 0 ){
	  qpPreStepEnergy_meV = theStep1.preStepKinEnergy_eV * 1000;
	  qpPostStepEnergy_meV = theStep1.postStepKinEnergy_eV * 1000;
	  isBogoliubovQP = true;
	}
	//If we're on the last step, take the final track point
	if( isBogoliubovQP && iS == theTrack.stepVect.size()-1 ){
	  qpLX = theStep1.postStepX_mm;
	  qpLY = theStep1.postStepY_mm;
	  qpLZ = theStep1.postStepZ_mm;
	}
      }
    }

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	int trackID = theStep1.trackID;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double e0 = theStep1.preStepKinEnergy_eV;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double e1 = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double dR = TMath::Power(x1*x1 + y1*y1 + z1*z1, 0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double stepDeltaT_ns = t1-t0;
	double dz_1 = z1-z0;
	std::string particleName = theStep1.particleName;
	std::string processName = theStep1.processName;

	//If we're looking at a QP, plot generic info
	if( particleName == "BogoliubovQP" ){
	  h_QP_postStepPointXY->Fill(x1,y1);
	  h_QP_postStepPointXZ->Fill(x1,z1);
	  h_QP_postStepPointYZ->Fill(y1,z1);
	  h_QP_deltaT->Fill(t1-t0);
	  h_QP_preStepEnergy->Fill(e0*1000);
	  h_QP_postStepEnergy->Fill(e1*1000);

	  //If we're looking at a QP *and* the step is also a qprecombination step, then a different
	  //set of histograms gets filled
	  if( processName == "bogoliubovQPRecombination" ){
	    h_QPrec_postStepPointXY->Fill(x1,y1);
	    h_QPrec_postStepPointXZ->Fill(x1,z1);
	    h_QPrec_postStepPointYZ->Fill(y1,z1);
	    h_QPrec_DeltaTVsPreStepEnergy->Fill(e0*1000,t1-t0);
	    h_QPrec_DeltaTVsPreStepEnergy_LowE->Fill(e0*1000,t1-t0);
	    h_QPrec_deltaT->Fill(t1-t0);
	    h_QPrec_preStepEnergy->Fill(e0*1000);
	    h_QPrec_postStepEnergy->Fill(e1*1000);
	  }
	  
	}

	//If we're looking at phonons...
	if( particleName == "phononTF" || particleName == "phononTS" || particleName == "phononL" ){

	  //Associate the polarization to a number
	  double polariz = -1;
	  if( particleName == "phononL" ) polariz = 0.5;
	  if( particleName == "phononTF" ) polariz = 1.5;
	  if( particleName == "phononTS" ) polariz = 2.5;
	  h_phononPolarization->Fill(polariz);	  
	  h_phononThetaPhi->Fill(TMath::ACos(dZ/dR)*180./TMath::Pi(),TMath::ATan2(dY,dX)*180./TMath::Pi());	  
	  
	  //And, if we're looking at specifically the trackID = 2 particle (a phonon), then we can use the existing energy
	  //of the initial QP to get a distribution of which energies come from which other energies. Can't do this easily for
	  //later tracks because bookkeeping will be a nightmare
	  if( trackID == 2 ){
	    h_phononEnergyVsQPPreStepEnergy->Fill(qpPreStepEnergy_meV,e0*1000);
	  }
	}
      }
    }
  }


  //Let's do some fits
  
  //More will go here in a bit
  //Set up fitter to do lots of iterations
  ROOT::Math::MinimizerOptions temp;
  temp.SetMaxIterations(10000);
  
  //Do fits to vertical slices of the deltaTVsPreStepEnergy plot, and get the decay time
  double fitMin_ns = 0.1;
  //double fitMax_ns = 1000000; //For high energy
  double fitMax_ns = 10000; //For high energy
  TF1 * expFit = new TF1("expFit","[0]*TMath::Exp(-1*x/[1])",fitMin_ns,fitMax_ns);
  TGraphErrors * g_fitExponents = new TGraphErrors();
  int offset = 0; //For plotting
  for( int iBX = 1; iBX < h_QPrec_DeltaTVsPreStepEnergy_LowE->GetNbinsX(); ++iBX){
    std::cout << " pre-step energy bin center: " << h_QPrec_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX) << std::endl;
    if( h_QPrec_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX) < 1*gap_meV ){
      std::cout << "skipping this round." << std::endl;
      offset += 1;
      continue;
    }
    
    char name[400];
    sprintf(name,"theSlice_%d",iBX);
    TH1F * theSlice = (TH1F*)h_QPrec_DeltaTVsPreStepEnergy_LowE->ProjectionY(name,iBX,iBX,"");
    theSlice->Rebin(20);
    TCanvas * c1 = new TCanvas();
    theSlice->Draw("HIST");
    expFit->SetParameter(0,theSlice->GetBinContent(1)); //Hardcoded start point (based on number of entries)
    expFit->SetParameter(1,theSlice->GetMean()); //Hardcoded start point for high-energy
    theSlice->Fit(expFit,"","",fitMin_ns,fitMax_ns);
    double fitExponent_ns = expFit->GetParameter(1);
    double fitExponentError_ns = expFit->GetParError(1);
    expFit->Draw("Lsame");
    
    std::cout << "Fit exponent: " << fitExponent_ns << std::endl;
    
    //Convert the fit exponent to taus/tau0
    double tauRDivTau0 = fitExponent_ns / tau0_qp_ns;
    double tauRDivTau0_error = tauRDivTau0 * (fitExponentError_ns / fitExponent_ns);
    std::cout << "TauRDivTau0 = " << tauRDivTau0 << ", and plotting." << std::endl;
    
    g_fitExponents->SetPoint(iBX-offset-1,h_QPrec_DeltaTVsPreStepEnergy_LowE->GetXaxis()->GetBinCenter(iBX),tauRDivTau0);
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


//----------------------------------------------------------------------------------------
//This is an analysis of converted steps to make sure we understand that phonon radiation by
//quasiparticles is working as expected
void AnalyzeConvertedSteps_BogoliubovQPDiffusionProcessStudy(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In BogoliubovQP Diffusion Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;

  double gapE_meV = 0.175157;
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPDiffusionProcessStudy.root","RECREATE");  

  //Defining variables
  int nStepsE = 10;
  double minE_meV = 0.0;
  double maxE_meV = 0.6;
  double deltaE_meV = (maxE_meV-minE_meV)/nStepsE;
  int nStepsR = 200;
  double minR_mm = 0.0;
  double maxR_mm = 40.0;
  int nStepsX = 100;
  double minX_mm = -40.0;
  double maxX_mm = 40.0;
  int nStepsY = 100;
  double minY_mm = -40.0;
  double maxY_mm = 40.0;
  int nStepsZ = 200;
  double minZ_mm = -12.8;
  double maxZ_mm = 25.0;
  int nStepsDeltaTLong = 4000;
  double minDeltaTLong_ns = 0;
  double maxDeltaTLong_ns = 40000000; //500E4
  int nStepsT = 40;
  double minT_ns = 0;
  double maxT_ns = 0.25e8;
  int nStepsTLong = 400;
  double maxTLong_ns = 1e9;


  
  int nStepsDeltaT = 5000;
  double minDeltaT_ns = 0;
  double maxDeltaT_ns = 500;
  
  //Defining generic QP histograms (all steps)
  TH2F * h_QP_postStepPointXY = new TH2F("h_QP_postStepPointXY","QP PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QP_postStepPointXZ = new TH2F("h_QP_postStepPointXZ","QP PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QP_postStepPointYZ = new TH2F("h_QP_postStepPointYZ","QP PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_QP_deltaT = new TH1F("h_QP_deltaT","QP Step DeltaT; DeltaT [ns]; Steps",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH1F * h_QP_preStepEnergy = new TH1F("h_QP_preStepKinEnergy","QP Pre-step Kinetic energy; Kinetic Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_postStepEnergy = new TH1F("h_QP_postStepKinEnergy","QP Post-step Kinetic energy; Kinetic Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH2F * h_QP_deltaT_postStepR = new TH2F("h_QP_deltaT_postStepR","QP Post-step R vs. DeltaT; Delta T [ns]; Post-Step R [mm];",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns,1000,0,40);
					
				       
  
  //Setting up energy dependent diffusion histograms
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPRVsTime;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPXVsTime;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPYVsTime;

  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPRVsTimeLong;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPXVsTimeLong;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPYVsTimeLong;  
  for( int iSE = 0; iSE < nStepsE; ++iSE ){
    double energyLow = iSE*deltaE_meV + minE_meV;
    double energyHigh = energyLow+deltaE_meV;
    std::pair<double,double> energyBin(energyLow,energyHigh);
    char nameR[400];
    char nameX[400];
    char nameY[400];
    sprintf(nameR,"h_QP_R_vs_AbsPreStepTime_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameX,"h_QP_X_vs_AbsPreStepTime_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameY,"h_QP_Y_vs_AbsPreStepTime_Energy_%.4f_%.4f",energyLow,energyHigh);    
    TH2F * h_QP_R_vs_AbsPreStepTime = new TH2F(nameR,"QP Radial position vs. absolute pre-step time; R [mm]; Pre-step T [ns]; nQPs",nStepsR,minR_mm,maxR_mm,nStepsT,minT_ns,maxT_ns);
    TH2F * h_QP_X_vs_AbsPreStepTime = new TH2F(nameX,"QP X position vs. absolute pre-step time; X [mm]; Pre-step T [ns]; nQPs",nStepsX,minX_mm,maxX_mm,nStepsT,minT_ns,maxT_ns);
    TH2F * h_QP_Y_vs_AbsPreStepTime = new TH2F(nameY,"QP Y position vs. absolute pre-step time; Y [mm]; Pre-step T [ns]; nQPs",nStepsY,minY_mm,maxY_mm,nStepsT,minT_ns,maxT_ns);
    map_energyBin_QPRVsTime.emplace(energyBin,h_QP_R_vs_AbsPreStepTime);
    map_energyBin_QPXVsTime.emplace(energyBin,h_QP_X_vs_AbsPreStepTime);
    map_energyBin_QPYVsTime.emplace(energyBin,h_QP_Y_vs_AbsPreStepTime);


    

    sprintf(nameR,"h_QP_R_vs_AbsPreStepTimeLong_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameX,"h_QP_X_vs_AbsPreStepTimeLong_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameY,"h_QP_Y_vs_AbsPreStepTimeLong_Energy_%.4f_%.4f",energyLow,energyHigh);    
    TH2F * h_QP_R_vs_AbsPreStepTimeLong = new TH2F(nameR,"QP Radial position vs. absolute pre-step time; R [mm]; Pre-step T [ns]; nQPs",nStepsR,minR_mm,maxR_mm,nStepsTLong,minT_ns,maxTLong_ns);
    TH2F * h_QP_X_vs_AbsPreStepTimeLong = new TH2F(nameX,"QP X position vs. absolute pre-step time; X [mm]; Pre-step T [ns]; nQPs",nStepsX,minX_mm,maxX_mm,nStepsTLong,minT_ns,maxTLong_ns);
    TH2F * h_QP_Y_vs_AbsPreStepTimeLong = new TH2F(nameY,"QP Y position vs. absolute pre-step time; Y [mm]; Pre-step T [ns]; nQPs",nStepsY,minY_mm,maxY_mm,nStepsTLong,minT_ns,maxTLong_ns);    
    map_energyBin_QPRVsTimeLong.emplace(energyBin,h_QP_R_vs_AbsPreStepTimeLong);
    map_energyBin_QPXVsTimeLong.emplace(energyBin,h_QP_X_vs_AbsPreStepTimeLong);
    map_energyBin_QPYVsTimeLong.emplace(energyBin,h_QP_Y_vs_AbsPreStepTimeLong);
  }
													     
  //Okay, time to start looping over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    if( iE % 1000 ) std::cout << "Event Number " << iE << std::endl;

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	int trackID = theStep1.trackID;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double e0 = theStep1.preStepKinEnergy_eV;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double e1 = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r0 = TMath::Power(x0*x0 + y0*y0,0.5);
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double dR = TMath::Power(x1*x1 + y1*y1 + z1*z1, 0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double stepDeltaT_ns = t1-t0;
	double dz_1 = z1-z0;
	std::string particleName = theStep1.particleName;
	std::string processName = theStep1.processName;

	//Double-check to make sure we're only looking at QPs
	if( particleName == "BogoliubovQP" ){
	  h_QP_postStepPointXY->Fill(x1,y1);
	  h_QP_postStepPointXZ->Fill(x1,z1);
	  h_QP_postStepPointYZ->Fill(y1,z1);
	  h_QP_deltaT->Fill(t1-t0);
	  h_QP_deltaT_postStepR->Fill(t1-t0,r1);
	  h_QP_preStepEnergy->Fill(e0*1000);
	  h_QP_postStepEnergy->Fill(e1*1000);

	  //Now use the energy to figure out which diffusion-vs-time histogram to fill
	  for(std::map<std::pair<double,double>,TH2F*>::iterator it = map_energyBin_QPRVsTime.begin(); it != map_energyBin_QPRVsTime.end(); ++it ){
	    double minEnergy = it->first.first;
	    double maxEnergy = it->first.second;
	    if( e0*1000 < maxEnergy && e0*1000 >= minEnergy ){
	      map_energyBin_QPRVsTime[it->first]->Fill(r0,t0);
	      map_energyBin_QPXVsTime[it->first]->Fill(x0,t0);
	      map_energyBin_QPYVsTime[it->first]->Fill(y0,t0);


	      map_energyBin_QPRVsTimeLong[it->first]->Fill(r0,t0);
	      map_energyBin_QPXVsTimeLong[it->first]->Fill(x0,t0);
	      map_energyBin_QPYVsTimeLong[it->first]->Fill(y0,t0);	      
	    }
	  }
	}
      }
    }
  }

  //Now do some fits, looping over the energy bins
  TGraphErrors * g1_bin4 = new TGraphErrors();
  TGraphErrors * g1_bin5 = new TGraphErrors();
  TGraphErrors * g1_bin6 = new TGraphErrors();
  TGraphErrors * g1_bin7 = new TGraphErrors();
  TGraphErrors * g1_bin8 = new TGraphErrors();
  TGraphErrors * g1_avg = new TGraphErrors();
  TGraphErrors * g1fit_bin4 = new TGraphErrors();
  TGraphErrors * g1fit_bin5 = new TGraphErrors();
  TGraphErrors * g1fit_bin6 = new TGraphErrors();
  TGraphErrors * g1fit_bin7 = new TGraphErrors();
  TGraphErrors * g1fit_bin8 = new TGraphErrors();
  TF1 * myGauss = new TF1("myGauss","[0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))",minX_mm,maxX_mm);
  int counter = 0;
  std::vector<TH1F*> storageVect;
  for( std::map<std::pair<double,double>,TH2F*>::iterator it = map_energyBin_QPYVsTime.begin(); it != map_energyBin_QPYVsTime.end(); ++it ){
    double minEnergy = it->first.first;
    double maxEnergy = it->first.second;
    double avgEnergy = (minEnergy + maxEnergy)/2.0;

    if( it->second->GetEntries() == 0 ) continue;
    
    //Get bins 4, 5,6, 7, and 8 and extract info
    TH1F * timeBin4 = (TH1F*)it->second->ProjectionX("bin4",4,4);
    TH1F * timeBin5 = (TH1F*)it->second->ProjectionX("bin5",5,5);
    TH1F * timeBin6 = (TH1F*)it->second->ProjectionX("bin6",6,6);
    TH1F * timeBin7 = (TH1F*)it->second->ProjectionX("bin7",7,7);
    TH1F * timeBin8 = (TH1F*)it->second->ProjectionX("bin8",8,8);

    //Get the time of the slices
    double timeOfBin4 = it->second->GetYaxis()->GetBinCenter(4);
    double timeOfBin5 = it->second->GetYaxis()->GetBinCenter(5);
    double timeOfBin6 = it->second->GetYaxis()->GetBinCenter(6);
    double timeOfBin7 = it->second->GetYaxis()->GetBinCenter(7);
    double timeOfBin8 = it->second->GetYaxis()->GetBinCenter(8);

    //Get the sigma from each slice
    double simpleSigmaBin4 = timeBin4->GetStdDev();
    double simpleSigmaBin5 = timeBin5->GetStdDev();
    double simpleSigmaBin6 = timeBin6->GetStdDev();
    double simpleSigmaBin7 = timeBin7->GetStdDev();
    double simpleSigmaBin8 = timeBin8->GetStdDev();

    //Compute from this sigma a diffusion constant, D (here referred to as Dn, but it's really D)
    double simpleDnBin4 = TMath::Power(simpleSigmaBin4,2)/2.0/timeOfBin4;
    double simpleDnBin5 = TMath::Power(simpleSigmaBin5,2)/2.0/timeOfBin5;
    double simpleDnBin6 = TMath::Power(simpleSigmaBin6,2)/2.0/timeOfBin6;
    double simpleDnBin7 = TMath::Power(simpleSigmaBin7,2)/2.0/timeOfBin7;
    double simpleDnBin8 = TMath::Power(simpleSigmaBin8,2)/2.0/timeOfBin8;

    //Basic standard error on standard deviation
    double uncSimpleSigmaBin4 = simpleSigmaBin4/sqrt(2.0*(timeBin4->GetEntries()-1));
    double uncSimpleSigmaBin5 = simpleSigmaBin5/sqrt(2.0*(timeBin5->GetEntries()-1));
    double uncSimpleSigmaBin6 = simpleSigmaBin6/sqrt(2.0*(timeBin6->GetEntries()-1));
    double uncSimpleSigmaBin7 = simpleSigmaBin7/sqrt(2.0*(timeBin7->GetEntries()-1));
    double uncSimpleSigmaBin8 = simpleSigmaBin8/sqrt(2.0*(timeBin8->GetEntries()-1));

    //Propagating above value to uncertainty in simpleDn formed from the standard deviation. Since D has a squared sigma in it,
    //we end up needing something not super trivial
    double uncDnBin4 = simpleDnBin4 * 2 * uncSimpleSigmaBin4 / simpleSigmaBin4;
    double uncDnBin5 = simpleDnBin5 * 2 * uncSimpleSigmaBin5 / simpleSigmaBin5;
    double uncDnBin6 = simpleDnBin6 * 2 * uncSimpleSigmaBin6 / simpleSigmaBin6;
    double uncDnBin7 = simpleDnBin7 * 2 * uncSimpleSigmaBin7 / simpleSigmaBin7;
    double uncDnBin8 = simpleDnBin8 * 2 * uncSimpleSigmaBin8 / simpleSigmaBin8;

    //Compute the average D from the five bins, and compute its uncertainty using the uncertainties above.
    double averageDn = (simpleDnBin4 + simpleDnBin5 + simpleDnBin6 + simpleDnBin7 + simpleDnBin8)/5.0;
    double uncAverageDn = TMath::Sqrt(TMath::Power(uncDnBin4/5.0,2)+TMath::Power(uncDnBin5/5.0,2)+TMath::Power(uncDnBin6/5.0,2)+TMath::Power(uncDnBin7/5.0,2)+TMath::Power(uncDnBin8/5.0,2));
    
    std::cout << "For bin with average energy: " << avgEnergy << " timeOfBin4: " << timeOfBin4 << ", simpleSigmaBin4: " << simpleSigmaBin4 << ", simpleDnBin4: " << simpleDnBin4 << std::endl;

    //Write points and errors. For now, kill the error on the individual measurements for clarity, so we can see the error bar associated with
    //the averaged graph.
    std::cout << "For bin with average energy: " << avgEnergy << " and [average energy/gap energy] = " << avgEnergy / gapE_meV << ", bin 2 shows Dn: " << simpleDnBin4 << ", bin5 shows Dn: " << simpleDnBin5 << ", bin6 shows Dn: " << simpleDnBin6 << std::endl;
    g1_bin4->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin4);
    g1_bin5->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin5);
    g1_bin6->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin6);
    g1_bin7->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin7);
    g1_bin8->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin8);
    g1_bin4->SetPointError(counter,0,0);//uncDnBin4);    
    g1_bin5->SetPointError(counter,0,0);//uncDnBin5);
    g1_bin6->SetPointError(counter,0,0);//uncDnBin6);
    g1_bin7->SetPointError(counter,0,0);//uncDnBin7);
    g1_bin8->SetPointError(counter,0,0);//uncDnBin8);
    g1_avg->SetPoint(counter,avgEnergy/gapE_meV,averageDn);
    g1_avg->SetPointError(counter,0,uncAverageDn);




    
    //Prepare TF1 for fits
    myGauss->SetParameter(0,10);
    myGauss->SetParameter(1,0);
    myGauss->SetParameter(2,3);
    
    //Drawing
    TCanvas * c0 = new TCanvas();
    timeBin4->SetLineColor(kRed);
    timeBin4->Draw("HIST");
    timeBin4->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin4->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin4->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin4);
    
    TCanvas * c1 = new TCanvas();
    timeBin5->SetLineColor(kBlue);
    timeBin5->Draw("HIST");
    timeBin5->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin5->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin5->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin5);

    
    TCanvas * c2 = new TCanvas();
    timeBin6->SetLineColor(kViolet);
    timeBin6->Draw("HIST");
    timeBin6->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin6->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin6->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin6);
    
    
    TCanvas * c3 = new TCanvas();
    timeBin7->SetLineColor(kGreen);
    timeBin7->Draw("HIST");
    timeBin7->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin7->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin7->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin7);

    
    TCanvas * c4 = new TCanvas();
    timeBin8->SetLineColor(kOrange);
    timeBin8->Draw("HIST");
    timeBin8->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin8->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin8->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin8);

    

    storageVect.push_back(timeBin4);
    storageVect.push_back(timeBin5);
    storageVect.push_back(timeBin6);
    storageVect.push_back(timeBin7);
    storageVect.push_back(timeBin8);
    counter++;
  }

  g1_bin4->SetMarkerColor(kRed);
  g1_bin4->SetMarkerStyle(21);
  g1fit_bin4->SetMarkerColor(kRed);
  g1fit_bin4->SetMarkerStyle(22);
  g1_bin5->SetMarkerColor(kOrange);
  g1_bin5->SetMarkerStyle(21);
  g1fit_bin5->SetMarkerColor(kOrange);
  g1fit_bin5->SetMarkerStyle(22);  
  g1_bin6->SetMarkerColor(kGreen);
  g1_bin6->SetMarkerStyle(21);
  g1fit_bin6->SetMarkerColor(kGreen);
  g1fit_bin6->SetMarkerStyle(22);  
  g1_bin7->SetMarkerColor(kBlue);
  g1_bin7->SetMarkerStyle(21);
  g1fit_bin7->SetMarkerColor(kBlue);
  g1fit_bin7->SetMarkerStyle(22);    
  g1_bin8->SetMarkerColor(kViolet);
  g1_bin8->SetMarkerStyle(21);
  g1fit_bin8->SetMarkerColor(kViolet);
  g1fit_bin8->SetMarkerStyle(22);      
  g1_bin4->GetXaxis()->SetTitle("QP Energy / Gap Energy");
  g1_bin4->GetYaxis()->SetTitle("Dn [mm^2/ns]");
  g1_avg->SetMarkerColor(kBlack);
  g1_avg->SetMarkerStyle(21);
  TCanvas * c1 = new TCanvas();
  g1_bin4->Draw("AP");
  g1_bin5->Draw("Psame");
  g1_bin6->Draw("Psame");
  g1_bin7->Draw("Psame");
  g1_bin8->Draw("Psame");
  g1_avg->Draw("Psame");

  //  g1fit_bin4->Draw("Psame");
  //g1fit_bin5->Draw("Psame");
  //g1fit_bin6->Draw("Psame");
  //g1fit_bin7->Draw("Psame");
  //g1fit_bin8->Draw("Psame");

  //Calculate the theory curve
  TGraph * g1_theory = new TGraph();
  for( int iR = 0; iR < 150; iR++ ){
    
    double ratio = 1.0/(1.0+iR*0.02);
    double Dn = 6E-6 * sqrt((1-ratio*ratio));
    g1_theory->SetPoint(iR,1.0/ratio,Dn);    
  }
  g1_theory->SetLineColor(kBlack);
  g1_theory->Draw("Lsame");

  //Create a legend
  TLegend * l1 = new TLegend();
  l1->AddEntry(g1_bin4,"Time Bin 4, Analytic");
  l1->AddEntry(g1_bin5,"Time Bin 5, Analytic");
  l1->AddEntry(g1_bin6,"Time Bin 6, Analytic");
  l1->AddEntry(g1_bin7,"Time Bin 7, Analytic");
  l1->AddEntry(g1_bin8,"Time Bin 8, Analytic");
  l1->AddEntry(g1_avg,"Average Dn, Analytic");
  //l1->AddEntry(g1fit_bin4,"Time Bin 4, Fit");
  //l1->AddEntry(g1fit_bin5,"Time Bin 5, Fit");
  //l1->AddEntry(g1fit_bin6,"Time Bin 6, Fit");
  //l1->AddEntry(g1fit_bin7,"Time Bin 7, Fit");
  //l1->AddEntry(g1fit_bin8,"Time Bin 8, Fit");
  
  l1->AddEntry(g1_theory,"Theory/Input");
  l1->Draw("same");






  


  
  outFile->Write();
}     




//----------------------------------------------------------------------------------------
//This is for understanding the boundary-limited diffusion case
void AnalyzeConvertedSteps_BogoliubovQPDiffusionProcessStudy_BoundaryLimitedDiffusion(std::string stepFileName)
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile(stepFileName.c_str());
  std::cout << "Done filling events." << std::endl;
  std::cout << "In BogoliubovQP Diffusion Process Study." << std::endl;
  std::cout << "Size of event vect: " << eventVect.size() << std::endl;

  double gapE_meV = 0.175157;
  
  //Output file
  TFile * outFile = new TFile("SteppingStudies_BogoliubovQPDiffusionProcessStudy_BoundaryLimitedDiffusion.root","RECREATE");  

  //Defining variables
  int nStepsE = 10;
  double minE_meV = 0.0;
  double maxE_meV = 0.6;
  double deltaE_meV = (maxE_meV-minE_meV)/nStepsE;
  int nStepsX = 100;
  double minX_mm = -40.0;
  double maxX_mm = 40.0;
  int nStepsY = 100;
  double minY_mm = -40.0;
  double maxY_mm = 40.0;
  int nStepsZ = 200;
  double minZ_mm = -12.8;
  double maxZ_mm = 25.0;
  int nStepsDeltaTLong = 4000;
  double minDeltaTLong_ns = 0;
  double maxDeltaTLong_ns = 40000000; //500E4
  int nStepsT = 40;
  double minT_ns = 0;
  double maxT_ns = 0.25e8;
  int nStepsTLong = 400;
  double maxTLong_ns = 1e9;


  
  int nStepsDeltaT = 5000;
  double minDeltaT_ns = 0;
  double maxDeltaT_ns = 500;
  
  //Defining generic QP histograms (all steps)
  TH2F * h_QP_postStepPointXY = new TH2F("h_QP_postStepPointXY","QP PostStepPoint XY; X [mm]; Y [mm];",1000,-40,40,1000,-40,40);
  TH2F * h_QP_postStepPointXZ = new TH2F("h_QP_postStepPointXZ","QP PostStepPoint XZ; X [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH2F * h_QP_postStepPointYZ = new TH2F("h_QP_postStepPointYZ","QP PostStepPoint YZ; Y [mm]; Z [mm];",1000,-40,40,1000,-12.8,25.0);
  TH1F * h_QP_deltaT = new TH1F("h_QP_deltaT","QP Step DeltaT; DeltaT [ns]; Steps",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns);
  TH1F * h_QP_preStepEnergy = new TH1F("h_QP_preStepKinEnergy","QP Pre-step Kinetic energy; Kinetic Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH1F * h_QP_postStepEnergy = new TH1F("h_QP_postStepKinEnergy","QP Post-step Kinetic energy; Kinetic Energy [meV]; Steps",nStepsE*100,minE_meV,maxE_meV);
  TH2F * h_QP_deltaT_postStepX = new TH2F("h_QP_deltaT_postStepX","QP Post-step X vs. DeltaT; Delta T [ns]; Post-Step X [mm];",nStepsDeltaTLong,minDeltaTLong_ns,maxDeltaTLong_ns,1000,-35,40);
					
				       
  
  //Setting up energy dependent diffusion histograms
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPXVsTime;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPYVsTime;

  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPXVsTimeLong;
  std::map<std::pair<double,double>,TH2F*> map_energyBin_QPYVsTimeLong;  
  for( int iSE = 0; iSE < nStepsE; ++iSE ){
    double energyLow = iSE*deltaE_meV + minE_meV;
    double energyHigh = energyLow+deltaE_meV;
    std::pair<double,double> energyBin(energyLow,energyHigh);
    char nameR[400];
    char nameX[400];
    char nameY[400];
    sprintf(nameX,"h_QP_X_vs_AbsPreStepTime_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameY,"h_QP_Y_vs_AbsPreStepTime_Energy_%.4f_%.4f",energyLow,energyHigh);    
    TH2F * h_QP_X_vs_AbsPreStepTime = new TH2F(nameX,"QP X position vs. absolute pre-step time; X [mm]; Pre-step T [ns]; nQPs",nStepsX,minX_mm,maxX_mm,nStepsT,minT_ns,maxT_ns);
    TH2F * h_QP_Y_vs_AbsPreStepTime = new TH2F(nameY,"QP Y position vs. absolute pre-step time; Y [mm]; Pre-step T [ns]; nQPs",nStepsY,minY_mm,maxY_mm,nStepsT,minT_ns,maxT_ns);
    map_energyBin_QPXVsTime.emplace(energyBin,h_QP_X_vs_AbsPreStepTime);
    map_energyBin_QPYVsTime.emplace(energyBin,h_QP_Y_vs_AbsPreStepTime);


    

    sprintf(nameX,"h_QP_X_vs_AbsPreStepTimeLong_Energy_%.4f_%.4f",energyLow,energyHigh);
    sprintf(nameY,"h_QP_Y_vs_AbsPreStepTimeLong_Energy_%.4f_%.4f",energyLow,energyHigh);    
    TH2F * h_QP_X_vs_AbsPreStepTimeLong = new TH2F(nameX,"QP X position vs. absolute pre-step time; X [mm]; Pre-step T [ns]; nQPs",nStepsX,minX_mm,maxX_mm,nStepsTLong,minT_ns,maxTLong_ns);
    TH2F * h_QP_Y_vs_AbsPreStepTimeLong = new TH2F(nameY,"QP Y position vs. absolute pre-step time; Y [mm]; Pre-step T [ns]; nQPs",nStepsY,minY_mm,maxY_mm,nStepsTLong,minT_ns,maxTLong_ns);
    map_energyBin_QPXVsTimeLong.emplace(energyBin,h_QP_X_vs_AbsPreStepTimeLong);
    map_energyBin_QPYVsTimeLong.emplace(energyBin,h_QP_Y_vs_AbsPreStepTimeLong);
  }
													     
  //Okay, time to start looping over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];

    if( iE % 1000 ) std::cout << "Event Number " << iE << std::endl;

    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      int maxNBounces = 0;
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep1 = theTrack.stepVect[iS];
	
	//Compute quantities
	int trackID = theStep1.trackID;
	double x0 = theStep1.preStepX_mm;
	double y0 = theStep1.preStepY_mm;
	double z0 = theStep1.preStepZ_mm;
	double e0 = theStep1.preStepKinEnergy_eV;
	double t0 = theStep1.preStepT_ns;
	double x1 = theStep1.postStepX_mm;
	double y1 = theStep1.postStepY_mm;
	double z1 = theStep1.postStepZ_mm;
	double e1 = theStep1.postStepKinEnergy_eV;
	double t1 = theStep1.postStepT_ns;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double r0 = TMath::Power(x0*x0 + y0*y0,0.5);
	double r1 = TMath::Power(x1*x1 + y1*y1,0.5);
	double dR = TMath::Power(x1*x1 + y1*y1 + z1*z1, 0.5);
	double stepLength1_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	double stepDeltaT_ns = t1-t0;
	double dz_1 = z1-z0;
	std::string particleName = theStep1.particleName;
	std::string processName = theStep1.processName;

	//Double-check to make sure we're only looking at QPs
	if( particleName == "BogoliubovQP" ){
	  h_QP_postStepPointXY->Fill(x1,y1);
	  h_QP_postStepPointXZ->Fill(x1,z1);
	  h_QP_postStepPointYZ->Fill(y1,z1);
	  h_QP_deltaT->Fill(t1-t0);
	  h_QP_deltaT_postStepX->Fill(t1-t0,x1);
	  h_QP_preStepEnergy->Fill(e0*1000);
	  h_QP_postStepEnergy->Fill(e1*1000);

	  //Now use the energy to figure out which diffusion-vs-time histogram to fill
	  for(std::map<std::pair<double,double>,TH2F*>::iterator it = map_energyBin_QPXVsTime.begin(); it != map_energyBin_QPXVsTime.end(); ++it ){
	    double minEnergy = it->first.first;
	    double maxEnergy = it->first.second;
	    if( e0*1000 < maxEnergy && e0*1000 >= minEnergy ){
	      map_energyBin_QPXVsTime[it->first]->Fill(x0,t0);
	      map_energyBin_QPYVsTime[it->first]->Fill(y0,t0);
	      map_energyBin_QPXVsTimeLong[it->first]->Fill(x0,t0);
	      map_energyBin_QPYVsTimeLong[it->first]->Fill(y0,t0);	      
	    }
	  }
	}
      }
    }
  }

  /*
  //Now do some fits, looping over the energy bins
  TGraphErrors * g1_bin4 = new TGraphErrors();
  TGraphErrors * g1_bin5 = new TGraphErrors();
  TGraphErrors * g1_bin6 = new TGraphErrors();
  TGraphErrors * g1_bin7 = new TGraphErrors();
  TGraphErrors * g1_bin8 = new TGraphErrors();
  TGraphErrors * g1_avg = new TGraphErrors();
  TGraphErrors * g1fit_bin4 = new TGraphErrors();
  TGraphErrors * g1fit_bin5 = new TGraphErrors();
  TGraphErrors * g1fit_bin6 = new TGraphErrors();
  TGraphErrors * g1fit_bin7 = new TGraphErrors();
  TGraphErrors * g1fit_bin8 = new TGraphErrors();
  TF1 * myGauss = new TF1("myGauss","[0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))",minX_mm,maxX_mm);
  int counter = 0;
  std::vector<TH1F*> storageVect;
  for( std::map<std::pair<double,double>,TH2F*>::iterator it = map_energyBin_QPXVsTime.begin(); it != map_energyBin_QPXVsTime.end(); ++it ){
    double minEnergy = it->first.first;
    double maxEnergy = it->first.second;
    double avgEnergy = (minEnergy + maxEnergy)/2.0;

    if( it->second->GetEntries() == 0 ) continue;
    
    //Get bins 4, 5,6, 7, and 8 and extract info
    TH1F * timeBin4 = (TH1F*)it->second->ProjectionX("bin4",4,4);
    TH1F * timeBin5 = (TH1F*)it->second->ProjectionX("bin5",5,5);
    TH1F * timeBin6 = (TH1F*)it->second->ProjectionX("bin6",6,6);
    TH1F * timeBin7 = (TH1F*)it->second->ProjectionX("bin7",7,7);
    TH1F * timeBin8 = (TH1F*)it->second->ProjectionX("bin8",8,8);

    //Get the time of the slices
    double timeOfBin4 = it->second->GetYaxis()->GetBinCenter(4);
    double timeOfBin5 = it->second->GetYaxis()->GetBinCenter(5);
    double timeOfBin6 = it->second->GetYaxis()->GetBinCenter(6);
    double timeOfBin7 = it->second->GetYaxis()->GetBinCenter(7);
    double timeOfBin8 = it->second->GetYaxis()->GetBinCenter(8);

    //Get the sigma from each slice
    double simpleSigmaBin4 = timeBin4->GetStdDev();
    double simpleSigmaBin5 = timeBin5->GetStdDev();
    double simpleSigmaBin6 = timeBin6->GetStdDev();
    double simpleSigmaBin7 = timeBin7->GetStdDev();
    double simpleSigmaBin8 = timeBin8->GetStdDev();

    //Compute from this sigma a diffusion constant, D (here referred to as Dn, but it's really D)
    double simpleDnBin4 = TMath::Power(simpleSigmaBin4,2)/2.0/timeOfBin4;
    double simpleDnBin5 = TMath::Power(simpleSigmaBin5,2)/2.0/timeOfBin5;
    double simpleDnBin6 = TMath::Power(simpleSigmaBin6,2)/2.0/timeOfBin6;
    double simpleDnBin7 = TMath::Power(simpleSigmaBin7,2)/2.0/timeOfBin7;
    double simpleDnBin8 = TMath::Power(simpleSigmaBin8,2)/2.0/timeOfBin8;

    //Basic standard error on standard deviation
    double uncSimpleSigmaBin4 = simpleSigmaBin4/sqrt(2.0*(timeBin4->GetEntries()-1));
    double uncSimpleSigmaBin5 = simpleSigmaBin5/sqrt(2.0*(timeBin5->GetEntries()-1));
    double uncSimpleSigmaBin6 = simpleSigmaBin6/sqrt(2.0*(timeBin6->GetEntries()-1));
    double uncSimpleSigmaBin7 = simpleSigmaBin7/sqrt(2.0*(timeBin7->GetEntries()-1));
    double uncSimpleSigmaBin8 = simpleSigmaBin8/sqrt(2.0*(timeBin8->GetEntries()-1));

    //Propagating above value to uncertainty in simpleDn formed from the standard deviation. Since D has a squared sigma in it,
    //we end up needing something not super trivial
    double uncDnBin4 = simpleDnBin4 * 2 * uncSimpleSigmaBin4 / simpleSigmaBin4;
    double uncDnBin5 = simpleDnBin5 * 2 * uncSimpleSigmaBin5 / simpleSigmaBin5;
    double uncDnBin6 = simpleDnBin6 * 2 * uncSimpleSigmaBin6 / simpleSigmaBin6;
    double uncDnBin7 = simpleDnBin7 * 2 * uncSimpleSigmaBin7 / simpleSigmaBin7;
    double uncDnBin8 = simpleDnBin8 * 2 * uncSimpleSigmaBin8 / simpleSigmaBin8;

    //Compute the average D from the five bins, and compute its uncertainty using the uncertainties above.
    double averageDn = (simpleDnBin4 + simpleDnBin5 + simpleDnBin6 + simpleDnBin7 + simpleDnBin8)/5.0;
    double uncAverageDn = TMath::Sqrt(TMath::Power(uncDnBin4/5.0,2)+TMath::Power(uncDnBin5/5.0,2)+TMath::Power(uncDnBin6/5.0,2)+TMath::Power(uncDnBin7/5.0,2)+TMath::Power(uncDnBin8/5.0,2));
    
    std::cout << "For bin with average energy: " << avgEnergy << " timeOfBin4: " << timeOfBin4 << ", simpleSigmaBin4: " << simpleSigmaBin4 << ", simpleDnBin4: " << simpleDnBin4 << std::endl;

    //Write points and errors. For now, kill the error on the individual measurements for clarity, so we can see the error bar associated with
    //the averaged graph.
    std::cout << "For bin with average energy: " << avgEnergy << " and [average energy/gap energy] = " << avgEnergy / gapE_meV << ", bin 2 shows Dn: " << simpleDnBin4 << ", bin5 shows Dn: " << simpleDnBin5 << ", bin6 shows Dn: " << simpleDnBin6 << std::endl;
    g1_bin4->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin4);
    g1_bin5->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin5);
    g1_bin6->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin6);
    g1_bin7->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin7);
    g1_bin8->SetPoint(counter,avgEnergy/gapE_meV,simpleDnBin8);
    g1_bin4->SetPointError(counter,0,0);//uncDnBin4);    
    g1_bin5->SetPointError(counter,0,0);//uncDnBin5);
    g1_bin6->SetPointError(counter,0,0);//uncDnBin6);
    g1_bin7->SetPointError(counter,0,0);//uncDnBin7);
    g1_bin8->SetPointError(counter,0,0);//uncDnBin8);
    g1_avg->SetPoint(counter,avgEnergy/gapE_meV,averageDn);
    g1_avg->SetPointError(counter,0,uncAverageDn);




    
    //Prepare TF1 for fits
    myGauss->SetParameter(0,10);
    myGauss->SetParameter(1,0);
    myGauss->SetParameter(2,3);
    
    //Drawing
    TCanvas * c0 = new TCanvas();
    timeBin4->SetLineColor(kRed);
    timeBin4->Draw("HIST");
    timeBin4->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin4->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin4->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin4);
    
    TCanvas * c1 = new TCanvas();
    timeBin5->SetLineColor(kBlue);
    timeBin5->Draw("HIST");
    timeBin5->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin5->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin5->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin5);

    
    TCanvas * c2 = new TCanvas();
    timeBin6->SetLineColor(kViolet);
    timeBin6->Draw("HIST");
    timeBin6->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin6->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin6->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin6);
    
    
    TCanvas * c3 = new TCanvas();
    timeBin7->SetLineColor(kGreen);
    timeBin7->Draw("HIST");
    timeBin7->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin7->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin7->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin7);

    
    TCanvas * c4 = new TCanvas();
    timeBin8->SetLineColor(kOrange);
    timeBin8->Draw("HIST");
    timeBin8->Fit(myGauss);
    myGauss->Draw("Lsame");
    //g1fit_bin8->SetPoint(counter,avgEnergy/gapE_meV,myGauss->GetParameter(2));
    //g1fit_bin8->SetPointError(counter,0,TMath::Power(myGauss->GetParError(2),2)/2.0/timeOfBin8);

    

    storageVect.push_back(timeBin4);
    storageVect.push_back(timeBin5);
    storageVect.push_back(timeBin6);
    storageVect.push_back(timeBin7);
    storageVect.push_back(timeBin8);
    counter++;
  }

  g1_bin4->SetMarkerColor(kRed);
  g1_bin4->SetMarkerStyle(21);
  g1fit_bin4->SetMarkerColor(kRed);
  g1fit_bin4->SetMarkerStyle(22);
  g1_bin5->SetMarkerColor(kOrange);
  g1_bin5->SetMarkerStyle(21);
  g1fit_bin5->SetMarkerColor(kOrange);
  g1fit_bin5->SetMarkerStyle(22);  
  g1_bin6->SetMarkerColor(kGreen);
  g1_bin6->SetMarkerStyle(21);
  g1fit_bin6->SetMarkerColor(kGreen);
  g1fit_bin6->SetMarkerStyle(22);  
  g1_bin7->SetMarkerColor(kBlue);
  g1_bin7->SetMarkerStyle(21);
  g1fit_bin7->SetMarkerColor(kBlue);
  g1fit_bin7->SetMarkerStyle(22);    
  g1_bin8->SetMarkerColor(kViolet);
  g1_bin8->SetMarkerStyle(21);
  g1fit_bin8->SetMarkerColor(kViolet);
  g1fit_bin8->SetMarkerStyle(22);      
  g1_bin4->GetXaxis()->SetTitle("QP Energy / Gap Energy");
  g1_bin4->GetYaxis()->SetTitle("Dn [mm^2/ns]");
  g1_avg->SetMarkerColor(kBlack);
  g1_avg->SetMarkerStyle(21);
  TCanvas * c1 = new TCanvas();
  g1_bin4->Draw("AP");
  g1_bin5->Draw("Psame");
  g1_bin6->Draw("Psame");
  g1_bin7->Draw("Psame");
  g1_bin8->Draw("Psame");
  g1_avg->Draw("Psame");

  //  g1fit_bin4->Draw("Psame");
  //g1fit_bin5->Draw("Psame");
  //g1fit_bin6->Draw("Psame");
  //g1fit_bin7->Draw("Psame");
  //g1fit_bin8->Draw("Psame");

  //Calculate the theory curve
  TGraph * g1_theory = new TGraph();
  for( int iR = 0; iR < 150; iR++ ){
    
    double ratio = 1.0/(1.0+iR*0.02);
    double Dn = 6E-6 * sqrt((1-ratio*ratio));
    g1_theory->SetPoint(iR,1.0/ratio,Dn);    
  }
  g1_theory->SetLineColor(kBlack);
  g1_theory->Draw("Lsame");

  //Create a legend
  TLegend * l1 = new TLegend();
  l1->AddEntry(g1_bin4,"Time Bin 4, Analytic");
  l1->AddEntry(g1_bin5,"Time Bin 5, Analytic");
  l1->AddEntry(g1_bin6,"Time Bin 6, Analytic");
  l1->AddEntry(g1_bin7,"Time Bin 7, Analytic");
  l1->AddEntry(g1_bin8,"Time Bin 8, Analytic");
  l1->AddEntry(g1_avg,"Average Dn, Analytic");
  //l1->AddEntry(g1fit_bin4,"Time Bin 4, Fit");
  //l1->AddEntry(g1fit_bin5,"Time Bin 5, Fit");
  //l1->AddEntry(g1fit_bin6,"Time Bin 6, Fit");
  //l1->AddEntry(g1fit_bin7,"Time Bin 7, Fit");
  //l1->AddEntry(g1fit_bin8,"Time Bin 8, Fit");
  
  l1->AddEntry(g1_theory,"Theory/Input");
  l1->Draw("same");
  */





  


  
  outFile->Write();
}     













  






















//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// old shit -- delete soon



void AnalyzeConvertedStepsPhononReflectionStudy()
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile("/Users/ryanlinehan/QSC/Sims/Geant4/phonon-build/StepInformationFile.txt");
  std::cout << "Done filling events." << std::endl;

  //Output file
  TFile * outFile = new TFile("SteppingStudiesPhononReflectionOutput.root","RECREATE");
  
  //Define plotting histograms
  TH1F * h_primaryEnergy = new TH1F("h_primaryEnergy","Primary Phonon Energy; Log10(Energy [eV]); Primaries",100,-5,-1);
  TH1F * h_nStepsInTrack = new TH1F("h_nStepsInTrack","Number of Steps in the Track; NSteps; Tracks",205,0,205);
  TH2F * h_stepLastXY = new TH2F("h_stepFinalX_stepFinalY","Step Final X, Y; X [mm]; Y [mm];",4000,-4.1,4.1,4000,-4.1,4.1);
  TH2F * h_stepLastXZ = new TH2F("h_stepFinalX_stepFinalZ","Step Final X, Z; X [mm]; Z [mm];",400,-4.1,4.1,400,4.6,5.1);
  TH2F * h_stepLastYZ = new TH2F("h_stepFinalY_stepFinalZ","Step Final Y, Z; Y [mm]; Z [mm];",400,-4.1,4.1,400,4.6,5.1);
  TH2F * h_stepLastXYTop = new TH2F("h_stepFinalX_stepFinalY_Top","Step Final X (Top of Chip), Y; X [mm]; Y [mm];",4000,-4.1,4.1,4000,-4.1,4.1);
  TH2F * h_stepLastXYBottom = new TH2F("h_stepFinalX_stepFinalY_Bottom","Step Final X (Bottom of Chip), Y; X [mm]; Y [mm];",4000,-4.1,4.1,4000,-4.1,4.1);

  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];
      h_nStepsInTrack->Fill(theTrack.stepVect.size());
      
      //Check to see if this track has more than one step
      if( theTrack.stepVect.size() > 1 ){
	std::cout << "This track in event " << theTrack.eventID << " has more than one step, amazingly." << std::endl;
      }
      
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep = theTrack.stepVect[iS];

	//Compute quantities
	double x0 = theStep.preStepX_mm;
	double y0 = theStep.preStepY_mm;
	double z0 = theStep.preStepZ_mm;
	double x1 = theStep.postStepX_mm;
	double y1 = theStep.postStepY_mm;
	double z1 = theStep.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double stepLength_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	std::string processName = theStep.processName;

	//Select only the last step for tracks that have multiple steps (i.e. scatters)
	if( iS != theTrack.stepVect.size()-1 ){ continue; }
	
	h_stepLastXY->Fill(x1,y1);
	h_stepLastXZ->Fill(x1,z1);
	h_stepLastYZ->Fill(y1,z1);
	
	//Now, choose only steps that land on the bottom face
	if( z1 < 4.63 ){
	  h_stepLastXYBottom->Fill(x1,y1);
	}
	if( z1 > 4.99 ){
	  h_stepLastXYTop->Fill(x1,y1);
	}
      }
    }
  }

  outFile->Write();
  
  
}





void AnalyzeConvertedStepsPhononFocusingStudy()
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile("/Users/ryanlinehan/QSC/Sims/Geant4/phonon-build/StepInformationFile.txt");
  std::cout << "Done filling events." << std::endl;

  //Output file
  TFile * outFile = new TFile("SteppingStudiesPhononFocusingOutput.root","RECREATE");
  
  //Define plotting histograms
  TH1F * h_primaryEnergy = new TH1F("h_primaryEnergy","Primary Phonon Energy; Log10(Energy [eV]); Primaries",100,-5,-1);
  TH2F * h_stepLastXY = new TH2F("h_stepFinalX_stepFinalY","Step Final X, Y; X [mm]; Y [mm];",400,-4.1,4.1,400,-4.1,4.1);
  TH2F * h_stepLastXZ = new TH2F("h_stepFinalX_stepFinalZ","Step Final X, Z; X [mm]; Z [mm];",400,-4.1,4.1,400,4.6,5.1);
  TH2F * h_stepLastYZ = new TH2F("h_stepFinalY_stepFinalZ","Step Final Y, Z; Y [mm]; Z [mm];",400,-4.1,4.1,400,4.6,5.1);
  TH2F * h_stepLastXYTop = new TH2F("h_stepFinalX_stepFinalY_Top","Step Final X (Top of Chip), Y; X [mm]; Y [mm];",400,-4.1,4.1,400,-4.1,4.1);
  TH2F * h_stepLastXYBottom = new TH2F("h_stepFinalX_stepFinalY_Bottom","Step Final X (Bottom of Chip), Y; X [mm]; Y [mm];",400,-4.1,4.1,400,-4.1,4.1);

  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){
      Track theTrack = theEvent.trackVect[iT];

      //Check to see if this track has more than one step
      if( theTrack.stepVect.size() > 1 ){
	std::cout << "This track in event " << theTrack.eventID << " has more than one step, amazingly." << std::endl;
      }
      
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep = theTrack.stepVect[iS];

	//Compute quantities
	double x0 = theStep.preStepX_mm;
	double y0 = theStep.preStepY_mm;
	double z0 = theStep.preStepZ_mm;
	double x1 = theStep.postStepX_mm;
	double y1 = theStep.postStepY_mm;
	double z1 = theStep.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double stepLength_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	std::string processName = theStep.processName;

	//Select only the last step for tracks that have multiple steps (i.e. scatters)
	if( iS != theTrack.stepVect.size()-1 ){ continue; }
	
	h_stepLastXY->Fill(x1,y1);
	h_stepLastXZ->Fill(x1,z1);
	h_stepLastYZ->Fill(y1,z1);
	
	//Now, choose only steps that land on the bottom face
	if( z1 < 4.63 ){
	  h_stepLastXYBottom->Fill(x1,y1);
	}
	if( z1 > 4.99 ){
	  h_stepLastXYTop->Fill(x1,y1);
	}
      }
    }
  }

  outFile->Write();
  
  
}



void AnalyzeConvertedStepsPhononDownconversionStudy()
{
  std::vector<Event> eventVect = ReadInG4CMPStepTextFile("/Users/ryanlinehan/QSC/Sims/Geant4/phonon-build/StepInformationFile.txt");
  std::cout << "Done filling events." << std::endl;

  //Output file
  TFile * outFile = new TFile("SteppingStudiesPhononDownconversionOutput.root","RECREATE");
  
  //Define plotting histograms
  TH1F * h_primaryEnergy = new TH1F("h_primaryEnergy","Primary Phonon Energy; Log10(Energy [eV]); Primaries",100,-5,-1);

  
  TH2F * h_primaryX_primaryY = new TH2F("h_primaryX_primaryY","Primary Photon Generation Location; X [mm]; Y [mm];",100,-8,8,100,-8,8);
  TH2F * h_primaryY_primaryZ = new TH2F("h_primaryY_primaryZ","Primary Photon Generation Location; Y [mm]; Z [mm];",100,-8,8,100,4.6,5.2);
  TH2F * h_primaryX_primaryZ = new TH2F("h_primaryX_primaryZ","Primary Photon Generation Location; X [mm]; Z [mm];",100,-8,8,100,4.6,5.2);
  TH1F * h_logStepLength = new TH1F("h_logStepLength","Step length; Log10(Step Length [mm]); Steps",200,-7,2);
  TH1F * h_linStepLength = new TH1F("h_linStepLength","Step length; Step Length [mm]; Steps",2000,0,5);
  TH2F * h_preStepEnergy_stepLength = new TH2F("h_preStepEnergy_stepLength","Step length vs. Pre-Step Energy; Log(Pre-Step Energy [eV]); Log10(Step Length [mm]);",200,-3.5,-0.5,200,-7,2);
  TH2F * h_postStepEnergy_stepLength = new TH2F("h_postStepEnergy_stepLength","Step length vs. Post-Step Energy; Log(Post-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);
  TH2F * h_preStepEnergy_postStepEnergy = new TH2F("h_preStepEnergy_postStepEnergy","Pre-step energy vs. Post-step energy; Log(Pre-Step Energy [eV]); Log(Post-Step Energy [eV])",200,-4,2,200,-4,2); //Sanity check
  
  //Stratified by phonon polarization
  TH2F * h_preStepEnergy_stepLength_long = new TH2F("h_preStepEnergy_stepLength_long","Step length vs. Pre-Step Energy, Longitudinal phonons; Log(Pre-Step Energy [eV]); Log10(Step Length [mm]);",200,-3.5,-0.5,200,-7,2);
  TH2F * h_postStepEnergy_stepLength_long = new TH2F("h_postStepEnergy_stepLength_long","Step length vs. Post-Step Energy, Longitudinal phonons; Log(Post-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);
  TH2F * h_preStepEnergy_stepLength_long_downconvert = new TH2F("h_preStepEnergy_stepLength_long_downconvert","Step length vs. Pre-Step Energy, Longitudinal phonons; Log(Pre-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);
  TH2F * h_preStepEnergy_stepLength_long_scatter = new TH2F("h_preStepEnergy_stepLength_long_scatter","Step length vs. Pre-Step Energy, Longitudinal phonons; Log(Pre-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);

  TH2F * h_preStepEnergy_stepLength_TF = new TH2F("h_preStepEnergy_stepLength_TF","Step length vs. Pre-Step Energy, Transverse Fast phonons; Log(Pre-Step Energy [eV]); Log10(Step Length [mm]);",200,-3.5,-0.5,200,-7,2);
  TH2F * h_postStepEnergy_stepLength_TF = new TH2F("h_postStepEnergy_stepLength_TF","Step length vs. Post-Step Energy, Transverse Fast phonons; Log(Post-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);

  TH2F * h_preStepEnergy_stepLength_TS = new TH2F("h_preStepEnergy_stepLength_TS","Step length vs. Pre-Step Energy, Transverse Slow phonons; Log(Pre-Step Energy [eV]); Log10(Step Length [mm]);",200,-3.5,-0.5,200,-7,2);
  TH2F * h_postStepEnergy_stepLength_TS = new TH2F("h_postStepEnergy_stepLength_TS","Step length vs. Post-Step Energy, Transverse Slow phonons; Log(Post-Step Energy [eV]; Log10(Step Length [mm]))",200,-3.5,-0.5,200,-7,2);

  //Mode-mixing studies
  TH1F * h_modeMixStepsLog = new TH1F("h_modeMixStepsLog","# Mode Mixing Steps vs. Energy; Log10(Energy [eV]); # Steps",200,-3.5,-0.5);
  TH1F * h_allStepsLog = new TH1F("h_allStepsLog","# All Steps vs. Energy; Log10(Energy [eV]); # Steps",200,-3.5,-0.5);
  TH1F * h_modeMixStepsLin = new TH1F("h_modeMixStepsLin","# Mode Mixing Steps vs. Energy; Log10(Energy [eV]); # Steps",200,0,0.1);
  TH1F * h_allStepsLin = new TH1F("h_allStepsLin","# All Steps vs. Energy; Log10(Energy [eV]); # Steps",200,0,0.1);
  TGraphErrors * modeMixingProbLog = new TGraphErrors(h_modeMixStepsLog->GetXaxis()->GetNbins());
  TGraphErrors * modeMixingProbLin = new TGraphErrors(h_modeMixStepsLog->GetXaxis()->GetNbins());
  
  
  
    
  
  //Loop over events
  for( int iE = 0; iE < eventVect.size(); ++iE ){
    Event theEvent = eventVect[iE];
    
    //Loop over tracks within events (treating all as basically equal)
    for( int iT = 0; iT < theEvent.trackVect.size(); ++iT ){

      //fill only with primary info
      Track theTrack = theEvent.trackVect[iT];
      if( theTrack.trackID == 1 ){
	h_primaryEnergy->Fill(TMath::Log10(theTrack.stepVect[0].preStepEnergy_eV));
	h_primaryX_primaryY->Fill(theTrack.stepVect[0].preStepX_mm,theTrack.stepVect[0].preStepY_mm);
	h_primaryY_primaryZ->Fill(theTrack.stepVect[0].preStepY_mm,theTrack.stepVect[0].preStepZ_mm);
	h_primaryX_primaryZ->Fill(theTrack.stepVect[0].preStepX_mm,theTrack.stepVect[0].preStepZ_mm);
      }
      
      //Loop over steps within tracks
      for( int iS = 0; iS < theTrack.stepVect.size(); ++iS ){
	Step theStep = theTrack.stepVect[iS];

	//Compute quantities
	double x0 = theStep.preStepX_mm;
	double y0 = theStep.preStepY_mm;
	double z0 = theStep.preStepZ_mm;
	double x1 = theStep.postStepX_mm;
	double y1 = theStep.postStepY_mm;
	double z1 = theStep.postStepZ_mm;
	double dX = x1-x0;
	double dY = y1-y0;
	double dZ = z1-z0;
	double stepLength_mm = TMath::Power(dX*dX + dY*dY + dZ*dZ,0.5 );
	std::string processName = theStep.processName;

	//If the process is transportation, ignore this step - it's boundary-limited.
	if( processName == "Transportation" ) continue;
	
	//Fill histograms
	h_preStepEnergy_stepLength->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	h_postStepEnergy_stepLength->Fill(TMath::Log10(theStep.postStepEnergy_eV),TMath::Log10(stepLength_mm));
	h_preStepEnergy_postStepEnergy->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(theStep.postStepEnergy_eV));
	h_linStepLength->Fill(stepLength_mm);
	h_logStepLength->Fill(TMath::Log10(stepLength_mm));

	//Phonon mode mixing - requires knowledge of the NEXT step as well - set a protection here, and need to go
	//two tracks back because the last step of many tracks is a Transport step, where the particle type doesn't
	//change by default (so the second-to-last step will "see" the Transportation step as the "iS+1" step)
	if( iS == theTrack.stepVect.size()-1 || iS == theTrack.stepVect.size()-2 ) continue; 
	if( theTrack.stepVect[iS].particleName != theTrack.stepVect[iS+1].particleName ){
	  h_modeMixStepsLin->Fill(theStep.preStepEnergy_eV);
	  h_modeMixStepsLog->Fill(TMath::Log10(theStep.preStepEnergy_eV));
	}
	h_allStepsLin->Fill(theStep.preStepEnergy_eV);
	h_allStepsLog->Fill(TMath::Log10(theStep.preStepEnergy_eV));
	

	
	
	if( theStep.particleName == "phononL" ){
	  h_preStepEnergy_stepLength_long->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	  h_postStepEnergy_stepLength_long->Fill(TMath::Log10(theStep.postStepEnergy_eV),TMath::Log10(stepLength_mm));
	  if( processName == "phononDownconversion" ){
	    h_preStepEnergy_stepLength_long_downconvert->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	    //	    std::cout << "downconvert log energy: " << TMath::Log10(theStep.postStepEnergy_eV) << ", log step length: " << TMath::Log10(stepLength_mm) << std::endl;
	  }
	  if( processName == "phononScattering" ){
	    h_preStepEnergy_stepLength_long_scatter->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	  }
	}
	if( theStep.particleName == "phononTF" ){
	  h_preStepEnergy_stepLength_TF->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	  h_postStepEnergy_stepLength_TF->Fill(TMath::Log10(theStep.postStepEnergy_eV),TMath::Log10(stepLength_mm));
	}
	if( theStep.particleName == "phononTS" ){
	  h_preStepEnergy_stepLength_TS->Fill(TMath::Log10(theStep.preStepEnergy_eV),TMath::Log10(stepLength_mm));
	  h_postStepEnergy_stepLength_TS->Fill(TMath::Log10(theStep.postStepEnergy_eV),TMath::Log10(stepLength_mm));
	}	
      }
    }
  }

  //Make the graphs showing probability of mode mixing
  for( int iB = 1; iB <= h_modeMixStepsLog->GetXaxis()->GetNbins(); ++iB ){
    double x = h_modeMixStepsLog->GetBinCenter(iB);
    double mixSteps = h_modeMixStepsLog->GetBinContent(iB);
    double allSteps = h_allStepsLog->GetBinContent(iB);
    double y = mixSteps/allSteps;
    if( allSteps == 0 ) y = 0; //To account for infinities
    double yErr = TMath::Power(mixSteps,0.5)/allSteps;

    modeMixingProbLog->SetPoint(iB-1,x,y);
    modeMixingProbLog->SetPointError(iB-1,0,yErr);
  }
  modeMixingProbLog->GetXaxis()->SetTitle("Log10(Phonon Energy [eV])");
  modeMixingProbLog->GetYaxis()->SetTitle("Probability of mode mixing");
  modeMixingProbLog->Write();
  
  //Make the graphs showing probability of mode mixing
  for( int iB = 1; iB <= h_modeMixStepsLin->GetXaxis()->GetNbins(); ++iB ){
    double x = h_modeMixStepsLin->GetBinCenter(iB);
    double mixSteps = h_modeMixStepsLin->GetBinContent(iB);
    double allSteps = h_allStepsLin->GetBinContent(iB);
    double y = mixSteps/allSteps;
    if( allSteps == 0 ) y = 0; //To account for infinities
    double yErr = TMath::Power(mixSteps,0.5)/allSteps;
    modeMixingProbLin->SetPoint(iB-1,x,y);
    modeMixingProbLin->SetPointError(iB-1,0,yErr);
  }
  modeMixingProbLin->GetXaxis()->SetTitle("Phonon Energy [eV]");
  modeMixingProbLin->GetYaxis()->SetTitle("Probability of mode mixing");
  modeMixingProbLin->Write();
  

  outFile->Write();
  
  /*
  std::cout << "Event 1:" << std::endl;
  for( int iT = 0; iT < eventVect[0].trackVect.size(); ++iT ){
    for( int iS = 0; iS < eventVect[0].trackVect[iT].stepVect.size(); ++iS ){
      std::cout << "Track: " << eventVect[0].trackVect[iT].trackID << ", Step. Pre-step X: " << eventVect[0].trackVect[iT].stepVect[iS].preStepX_mm
		<< ", pre-step Y: " << eventVect[0].trackVect[iT].stepVect[iS].preStepY_mm
		<< ", post-step X: " << eventVect[0].trackVect[iT].stepVect[iS].postStepY_mm << std::endl;
    }
  }
  */
  
}




//---------------------------------------------------------------------------------------
// Parsing function
std::vector<Event> ReadInG4CMPStepTextFile(std::string filename)
{
  std::vector<Event> output;
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  int eventID = -1;
  int runID = -1;
  int trackID = -1;
  int counter = 0;
  while(1){
    if(!infile.good()) break;
    if(infile.is_open()){
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between commas)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while(getline(check1,token,' ')){
	tokens.push_back(token);
      }
      if( tokens.size() == 0 ) break;


      //If our event or run number changes, then start a new event and push it back into the output
      bool startNewTrack = false;
      if( std::atoi(tokens[0].c_str()) != runID || std::atoi(tokens[1].c_str()) != eventID ){
	Event newEvent;
	newEvent.runID = std::atoi(tokens[0].c_str());
	newEvent.eventID = std::atoi(tokens[1].c_str());
	output.push_back(newEvent);
	runID = newEvent.runID;
	eventID = newEvent.eventID;
	startNewTrack = true;
	
	if( eventID % 10 == 0 ) std::cout << "Done reading " << eventID << " events." << std::endl;
	counter++;
      }

      
      //If our track number changes (or if we roll over to a new event), then start a new track and push back into the current event
      if( std::atoi(tokens[2].c_str()) != trackID || startNewTrack ){
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
