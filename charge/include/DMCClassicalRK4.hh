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
//
// class G4ClassicalRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using the classical 4th Runge-Kutta method.

// History:
// - Created: J.Apostolakis, V.Grichine - 30.1.97
// - Moved into G4MagErrorStepper: W.Wander <wwc@mit.edu> - 12/09/97
// -------------------------------------------------------------------

#include "G4AffineTransform.hh"
#include "G4MagErrorStepper.hh"
#include "EqEMFieldXtal.hh"

class DMCClassicalRK4 : public G4MagErrorStepper 
{
  public:  // with description
  G4AffineTransform normalToValley;
  G4AffineTransform valleyToNormal;

    DMCClassicalRK4(EqEMFieldXtal *EquationMotion, G4int numberOfVariables 
= 6) ;

    ~DMCClassicalRK4() ;

    // A stepper that does not know about errors.
    // It is used by the MagErrorStepper stepper.
   
    void DumbStepper( const G4double  yIn[],
                      const G4double  dydx[],
                            G4double  h,
                            G4double  yOut[] ) ;
      // Given values for the variables y[0,..,n-1] and their derivatives
      // dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
      // method to advance the solution over an interval h and return the
      // incremented variables as yout[0,...,n-1], which not be a distinct
      // array from y. The user supplies the routine RightHandSide(x,y,dydx),
      // which returns derivatives dydx at x. The source is routine rk4 from
      // NRC p. 712-713 .

  public:  // without description

    G4int IntegratorOrder() const { return 4; }
    void SetValleyTransform(G4AffineTransform xform);

  private:

    void StepWithEst( const G4double  yIn[],
                      const G4double  dydx[],
                            G4double  h,
                            G4double  yOut[],
                            G4double& alpha2,
                            G4double& beta2,
                      const G4double B1[],
                            G4double B2[] );
      // No longer used. Obsolete.

    DMCClassicalRK4(const DMCClassicalRK4&);
    DMCClassicalRK4& operator=(const DMCClassicalRK4&);
      // Private copy constructor and assignment operator.
  private:

    // G4int fNumberOfVariables ; // is set default to 6 in constructor

    G4double *dydxm, *dydxt, *yt; // scratch space - not state 
};
