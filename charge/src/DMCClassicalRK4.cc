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
//
// $Id$
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "DMCClassicalRK4.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"

//////////////////////////////////////////////////////////////////
//
// Constructor sets the number of variables (default = 6)

DMCClassicalRK4::DMCClassicalRK4(EqEMFieldXtal* EqRhs, G4int numberOfVariables) 
: G4MagErrorStepper(EqRhs, numberOfVariables)
{
   unsigned int noVariables= std::max(numberOfVariables,8); // For Time .. 7+1
 
   dydxm = new G4double[noVariables];
   dydxt = new G4double[noVariables]; 
   yt    = new G4double[noVariables]; 
   //normalToValley = EqRhs->GetNormalToValleyTransform();
   //valleyToNormal = EqRhs->GetValleyToNormalTransform();
   normalToValley = EqRhs->GetValleyTransform();
   valleyToNormal = normalToValley.Inverse();
}

////////////////////////////////////////////////////////////////
//
// Destructor

DMCClassicalRK4::~DMCClassicalRK4()
{
  delete[] dydxm;
  delete[] dydxt;
  delete[] yt;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for the variables y[0,..,n-1] and their derivatives
// dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
// method to advance the solution over an interval h and return the
// incremented variables as yout[0,...,n-1], which not be a distinct
// array from y. The user supplies the routine RightHandSide(x,y,dydx),
// which returns derivatives dydx at x. The source is routine rk4 from
// NRC p. 712-713 .

void
DMCClassicalRK4::DumbStepper( const G4double  yIn[],
                             const G4double  dydx[],
                                   G4double  h,
                                   G4double  yOut[])
{
   // G4cout << "DumbStepper" << G4endl;
  const G4int nvar = this->GetNumberOfVariables();   //  fNumberOfVariables(); 
  G4int i;

  // Initialise time to t0, needed when it is not updated by the integration.
  //        [ Note: Only for time dependent fields (usually electric) 
  //                  is it neccessary to integrate the time.] 
  G4ThreeVector xIn = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  //G4cout << "xIn: "<< xIn/cm << G4endl;
  G4ThreeVector pIn = G4ThreeVector(yIn[3], yIn[4], yIn[5]);
  
  G4ThreeVector xt;
  G4ThreeVector pt;

  for (i=0;i<3;i++)
  {
    xt[i] = xIn[i] + h/2.0*dydx[i];
    pt[i] = pIn[i] + h/2.0*dydx[i+3];
  }
  
  for(i=0;i<3;i++)
  {
    yt[i] = xt[i];
    yt[i+3] = pt[i];
  }
  RightHandSide(yt,dydxt) ;                   // 2nd Step K2=h*dydxt
  
  for (i=0;i<3;i++)
  {
    xt[i] = xIn[i] + h/2.0*dydxt[i];
    pt[i] = pIn[i] + h/2.0*dydxt[i+3];
  }
  
  for(i=0;i<3;i++)
  {
    yt[i] = xt[i];
    yt[i+3] = pt[i];
  }
  RightHandSide(yt,dydxm) ;                   // 3rd Step K3=h*dydxm
  
  for (i=0;i<3;i++)
  {
    xt[i] = xIn[i] + h/2.0*dydxm[i];
    pt[i] = pIn[i] + h/2.0*dydxm[i+3];
    dydxm[i] += dydxt[i];
    dydxm[i+3] += dydxt[i+3];
  }
  dydxm[7] += dydxt[7];
  
  for(i=0;i<3;i++)
  {
    yt[i] = xt[i];
    yt[i+3] = pt[i];
  }
  RightHandSide(yt,dydxt) ;                   // 4th Step K4=h*dydxt
 
  for(i=0;i<3;i++)    // Final RK4 output
  {
    xt[i] = xIn[i] + h/6.0*(dydx[i] + dydxt[i]+2.0*dydxm[i]);
    pt[i] = pIn[i] + h/6.0*(dydx[i+3] + dydxt[i+3]+2.0*dydxm[i+3]);
  }
  
  for (i=0;i<3;i++)
  {
    yOut[i] = xt[i];
    yOut[i+3] = pt[i];
  }
  yOut[7] = yIn[7] + h/6.0*(dydx[7]+dydxt[7]+2.0*dydxm[7]);
  if ( nvar == 12 )  { NormalisePolarizationVector ( yOut ); }
   /*   
    G4ThreeVector v_valley;
    v_valley[0] = pIn_HV[0]/.081/me;
    v_valley[1] = pIn_HV[1]/.081/me;
    v_valley[2] = pIn_HV[2]/1.58/me;
    G4ThreeVector v = valleyToNormal.TransformPoint(v_valley);
    v_valley = normalToValley.TransformPoint(pt);
    v_valley[0] /= .081*me;
    v_valley[1] /= .081*me;
    v_valley[2] /= 1.58*me;
    G4ThreeVector vOut = valleyToNormal.TransformPoint(v_valley);
    G4double t = h/6.0*(dydx[7]+dydxt[7]+2.0*dydxm[7]);
  G4cout << "Step Size: " << h/m << G4endl;
  G4cout << "vIn:       " << v/c_light*s/m << G4endl;
  G4cout << "vOut:      " << vOut/c_light*s/m << G4endl;
  G4cout << "t:         " << t/s << G4endl;
  */
  
}  // end of DumbStepper ....................................................

////////////////////////////////////////////////////////////////////
//
// StepWithEst

void
DMCClassicalRK4::StepWithEst( const G4double*,
                             const G4double*,
                                   G4double,
                                   G4double*,
                                   G4double&,
                                   G4double&,
                             const G4double*,
                                   G4double*  ) 
{
  G4Exception("DMCClassicalRK4::StepWithEst()", "GeomField0001",
              FatalException, "Method no longer used.");

}  // end of StepWithEst ......................................................

void DMCClassicalRK4::SetValleyTransform(G4AffineTransform xform){
  normalToValley = xform;
  valleyToNormal = xform.Inverse();
}
