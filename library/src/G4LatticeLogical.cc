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
/// \file exoticphysics/phonon/src/XLogicalLattice.cc
/// \brief Implementation of the XLogicalLattice class
//
// $Id$
//
#include "XLogicalLattice.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include <cmath>

XLogicalLattice::XLogicalLattice(){
  fVresTheta=0;
  fVresPhi=0;
  fDresTheta=0;
  fDresPhi=0;
  fA=0;
  fB=0;
  fDosL=0;
  fDosST=0;
  fDosFT=0;
  fBeta=0;
  fLambda=0;
  fGamma=0;
  fMu=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XLogicalLattice::~XLogicalLattice(){;}


void XLogicalLattice::SetDynamicalConstants(double Beta,double Gamma, double Lambda, double Mu)
{
  fBeta=Beta;
  fGamma=Gamma;
  fLambda=Lambda;
  fMu=Mu;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetScatteringConstant(G4double b)
{
  //Constant governing the rate of isotope scattering, for use with
  //XPhononScatteringProcess
  fB=b;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetAnhDecConstant(G4double a)
{
  //Constant governing rate of anharmonic down conversion of L-phonons,
  //for use with XPhononDownconversionProcess
  fA=a;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalLattice::SetLDOS(double LDOS)
{
  //Longitudinal phonon density of states
  fDosL=LDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetSTDOS(double STDOS)
{
  //Slow transverse phonon density of states
  fDosST=STDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetFTDOS(double FTDOS)
{
  //Fast transverse phonon density of states
  fDosFT=FTDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetBeta()
{
  return fBeta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetGamma()
{
  return fGamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetLambda()
{
  return fLambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetMu()
{
  return fMu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XLogicalLattice::GetScatteringConstant()
{
  return fB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XLogicalLattice::GetAnhDecConstant()
{
  return fA;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetLDOS()
{
  return fDosL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetSTDOS()
{
  return fDosST;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetFTDOS()
{
  return fDosFT;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


///////////////////////////////////////////
//Load map of group velocity scalars (m/s)
////////////////////////////////////////////
bool XLogicalLattice::LoadMap(int tRes, int pRes, int polarizationState, string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V fMap exceeds maximum resolution of "<<MAXRES<<" by "<<MAXRES<<". terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  fMapFile.clear();

  fThetaRes=tRes;
  fPhiRes=pRes;
  fMapFile.open(m.data());
  if(!fMapFile.is_open()) return false;
  
  for(int theta = 0; theta<fThetaRes; theta++){
    for(int phi = 0; phi<fPhiRes; phi++){
        fMapFile>>fMap[polarizationState][theta][phi];        
        //G4cout<<"\nXLogicalLattice::LoadMap:loading map "<<m<<" theta: "<<theta<<" phi: "<<phi<<" map: "<<fMap[polarizationState][theta][phi];
    }
  }
  fMapFile.close();
  G4cout<<"\nXLogicalLattice::LoadMap() sucessful (Vg scalars)."<<endl;
  fVresTheta=tRes; //store map dimensions
  fVresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


////////////////////////////////////
//Load map of group velocity unit vectors
///////////////////////////////////
bool XLogicalLattice::Load_NMap(int tRes, int pRes, int polarizationState, string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V map exceeds maximum resolution of "<<MAXRES<<" by "<<MAXRES<<" terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  fMapFile.clear();

  fThetaRes=tRes;
  fPhiRes=pRes;
  fMapFile.open(m.data());
  if(!fMapFile.is_open()) return false;
  
  for(int theta = 0; theta<fThetaRes; theta++){
    for(int phi = 0; phi<fPhiRes; phi++){
      for(int coord = 0; coord<3; coord++){
        fMapFile>>fN_map[polarizationState][theta][phi][coord];
      }
      //G4cout<<"\n theta: " <<theta<<" phi: "<<phi<<" [xyz]: "<<fN_map[polarizationState][theta][phi][0]<<" "<<fN_map[polarizationState][theta][phi][1]<<" "<<fN_map[polarizationState][theta][phi][2];
    }
  }
  G4cout<<"\n";
  fMapFile.close();
  G4cout<<"\nXLogicalLattice::Load_NMap() sucessful"<<endl;
  fDresTheta=tRes; //store map dimesnions
  fDresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::MapKtoV(int polarizationState,G4ThreeVector k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon velocity in m/s
    
  int useExtrapolation = 1; //uses bilinear extrapolation if ==1

  double theta, phi, tRes, pRes;
    
  int floorIntTheta, floorIntPhi, ceilIntTheta, ceilIntPhi, floorIntPhi2, ceilIntPhi2; //indices required for linear extrapolation
  double thetaFloor, thetaCeil, phiFloor, phiCeil, vFF, vCF, vFC, vCC; //values of theta, phi, v corresponding to floor/ceil integers

  tRes=pi/(fVresTheta);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fVresPhi);
  
  theta=k.getTheta();
  phi=k.getPhi();

  if(phi<0) phi = phi + 2*pi;
  if(theta>pi) theta=theta-pi;
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;

  if(useExtrapolation==1){
        floorIntTheta = (int) floor(theta/tRes);
        floorIntPhi = (int) floor(phi/pRes);
        ceilIntTheta = (int) ceil(theta/tRes);
        ceilIntPhi = (int) ceil(phi/pRes);
        
        if(ceilIntTheta==0) ceilIntTheta++; 
        if(ceilIntPhi==0) ceilIntPhi++; 
        
        thetaFloor = floorIntTheta*tRes;
        thetaCeil = ceilIntTheta*tRes;
        phiFloor = floorIntPhi*pRes;
        phiCeil = ceilIntPhi*pRes; 
        
        //velocity is zero when phi or theta reaches maximum map index (fDresTheta, fDresPhi)
        
        //if phi is greater than maximum phi in map, take coordinates (thetaCeil, 2*pi) as ceiling
        if(phi>(fDresPhi-1)*pRes){ 
            ceilIntPhi=0;
            phiCeil=2*M_PI;
        }
        
        //if theta is greater than maximum theta in map, take coordiates (theta, phi+PI) as upper bound
        if(theta>(fDresTheta-1)*tRes){ //(theta>dresTheta*tRes)
            
            ceilIntTheta=floorIntTheta;
            thetaCeil = thetaFloor;
            
            if(phi<M_PI){
                floorIntPhi2 = floor(phi+M_PI);
                ceilIntPhi2 = ceil(phi+M_PI);
            }
            else{
                floorIntPhi2 = floor(phi-M_PI);
                ceilIntPhi2 = ceil(phi-M_PI);
            }
            
            // perform bilinear interpolation
            vFF = fMap[polarizationState][floorIntTheta][floorIntPhi];
            vFC = fMap[polarizationState][floorIntTheta][ceilIntPhi];
            vCF = fMap[polarizationState][ceilIntTheta][floorIntPhi2];
            vCC = fMap[polarizationState][ceilIntTheta][ceilIntPhi2];
            
            if(vFF==0 || vFC==0 || vCF==0 || vCC==0){
                G4cout<<"\nvff vfc vcf vcc = "<<vFF<<", "<<vFC<<", "<<vCF<<", "<<vCC<<endl;
                G4cout<<"\n theta: " <<theta<<" phi: "<<phi<<endl;
            }
            
            double v_interp = (vFF/(2*pRes*tRes))*(M_PI+tRes-theta)*(phiCeil-phi) + (vCF/(2*pRes*tRes))*(theta-thetaFloor)*(phiCeil-phi)
            + (vFC/(2*pRes*tRes))*(thetaCeil-theta)*(phi-phiFloor) + (vCC/(2*pRes*tRes))*(theta-thetaFloor)*(phi-phiFloor);
            
            return v_interp;            
            
        }
        
        else{     
            // perform bilinear interpolation as normal
            vFF = fMap[polarizationState][floorIntTheta][floorIntPhi];
            vFC = fMap[polarizationState][floorIntTheta][ceilIntPhi];
            vCF = fMap[polarizationState][ceilIntTheta][floorIntPhi];
            vCC = fMap[polarizationState][ceilIntTheta][ceilIntPhi];
            
            
            double v_interp = (vFF/(pRes*tRes))*(thetaCeil-theta)*(phiCeil-phi) + (vCF/(pRes*tRes))*(theta-thetaFloor)*(phiCeil-phi)
            + (vFC/(pRes*tRes))*(thetaCeil-theta)*(phi-phiFloor) + (vCC/(pRes*tRes))*(theta-thetaFloor)*(phi-phiFloor);
            
            return v_interp;
        }
        
    }//end if(useExtrapolation==1)
    
    
    else{//use unmodified code if useExtrapolation not equal to 1
        
        if(int(theta/tRes)==fDresTheta) return fMap[polarizationState][int(theta/tRes)-1][int(phi/pRes)];
        
        else return fMap[polarizationState][int(theta/tRes)][int(phi/pRes)];
    }//end else
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4ThreeVector XLogicalLattice::MapKtoVDir(int polarizationState,G4ThreeVector k)
{  
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon propagation direction as dimensionless unit vector

  int useExtrapolation = 1; //uses linear extrapolation if ==1
    
  double theta, phi, tRes, pRes;
  

  tRes=pi/(fDresTheta-1);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fDresPhi-1);

  theta=k.getTheta();
  phi=k.getPhi(); 


  if(theta>pi) theta=theta-pi;
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;
  if(phi<0) phi = phi + 2*pi;

    if(useExtrapolation==1){
        double direction0 = BilinearInterpolateVDir(theta,phi,tRes,pRes,polarizationState,0);
        double direction1 = BilinearInterpolateVDir(theta,phi,tRes,pRes,polarizationState,1);
        double direction2 = BilinearInterpolateVDir(theta,phi,tRes,pRes,polarizationState,2);
        
        G4ThreeVector v(direction0,direction1,direction2);
        
        return v.unit();
    }
    
    else{//use code without interpolation ifunless useExtrapolation not == 1
        
        G4int iTheta = int(theta/tRes+0.5);
        G4int iPhi = int(phi/pRes+0.5);
        
        G4ThreeVector v(fN_map[polarizationState][iTheta][iPhi][0],
                        fN_map[polarizationState][iTheta][iPhi][1],
                        fN_map[polarizationState][iTheta][iPhi][2]);
        
        return v.unit();
                
        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



double XLogicalLattice::BilinearInterpolateVDir(double theta,double phi,double tRes,double pRes,int polarizationState,int mapIndex)
{
    //given the direction of the wave vector and the polarization state, looks up the velocity component for each dimension (mapIndex = 0,1,2)
    //and performs bilinear interpolation between values of the lookup table
    //returns the component of the velocity vector for input dimension (mapIndex = 0,1,2)
    
    int floorIntTheta, floorIntPhi, ceilIntTheta, ceilIntPhi, floorIntPhi2, ceilIntPhi2;
    double thetaFloor, thetaCeil, phiFloor, phiCeil, vDirFF, vDirFC, vDirCF, vDirCC;
    
    floorIntTheta = (int) floor(theta/tRes);
    floorIntPhi = (int) floor(phi/pRes);
    ceilIntTheta = (int) ceil(theta/tRes);
    ceilIntPhi = (int) ceil(phi/pRes);
    
    if(ceilIntTheta==0) ceilIntTheta++;
    if(ceilIntPhi==0) ceilIntPhi++; 
    
    thetaFloor = floorIntTheta*tRes;
    thetaCeil = ceilIntTheta*tRes;
    phiFloor = floorIntPhi*pRes;
    phiCeil = ceilIntPhi*pRes; 
    
    //velocity is zero when phi or theta reaches maximum map index (fDresTheta, fDresPhi)
    
    //if phi is greater than maximum phi in map, take coordinates (thetaCeil, 2*pi) as ceiling
    if(phi>(fDresPhi-1)*pRes){ 
        ceilIntPhi=0;
        phiCeil=2*M_PI;
    }
    
    
    //if theta is greater than maximum theta in map, take coordiates (theta, phi+PI) as upper bound 
    if(theta>(fDresTheta-1)*tRes){ 
        
        ceilIntTheta=floorIntTheta;
        thetaCeil = thetaFloor;
        
        if(phi<M_PI){
            floorIntPhi2 = floor(phi+M_PI);
            ceilIntPhi2 = ceil(phi+M_PI);
        }
        else{
            floorIntPhi2 = floor(phi-M_PI);
            ceilIntPhi2 = ceil(phi-M_PI);
        }
        
        // perform bilinear interpolation
        vDirFF = fN_map[polarizationState][floorIntTheta][floorIntPhi][mapIndex];
        vDirFC = fN_map[polarizationState][floorIntTheta][ceilIntPhi][mapIndex];
        vDirCF = fN_map[polarizationState][ceilIntTheta][floorIntPhi2][mapIndex];
        vDirCC = fN_map[polarizationState][ceilIntTheta][ceilIntPhi2][mapIndex];
        
        double vDir_interp = (vDirFF/(2*pRes*tRes))*(M_PI+tRes-theta)*(phiCeil-phi) + (vDirCF/(2*pRes*tRes))*(theta-thetaFloor)*(phiCeil-phi)
        + (vDirFC/(2*pRes*tRes))*(thetaCeil-theta)*(phi-phiFloor) + (vDirCC/(2*pRes*tRes))*(theta-thetaFloor)*(phi-phiFloor);
        
        return vDir_interp;  
        
    }
    
    else{ 
        
        // perform bilinear interpolation as usual
        vDirFF = fN_map[polarizationState][floorIntTheta][floorIntPhi][mapIndex];
        vDirFC = fN_map[polarizationState][floorIntTheta][ceilIntPhi][mapIndex];
        vDirCF = fN_map[polarizationState][ceilIntTheta][floorIntPhi][mapIndex];
        vDirCC = fN_map[polarizationState][ceilIntTheta][ceilIntPhi][mapIndex];
        
        double vDir_interp = (vDirFF/(pRes*tRes))*(thetaCeil-theta)*(phiCeil-phi) + (vDirCF/(pRes*tRes))*(theta-thetaFloor)*(phiCeil-phi)
        + (vDirFC/(pRes*tRes))*(thetaCeil-theta)*(phi-phiFloor) + (vDirCC/(pRes*tRes))*(theta-thetaFloor)*(phi-phiFloor);
        
        return vDir_interp; 
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

