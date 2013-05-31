//
// MODIFIED::
// $Id: G4EqMagElectricField.cc,v 1.14 2008/04/24 12:33:35 tnikitin Exp $
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
//  10.11.98   V.Grichine
//
// -------------------------------------------------------------------
// This version modified to allow oblique electron propagation
// D. Brandt
// 05 / 23 / 2011




#include "EqEMFieldXtal.hh"
#include "globals.hh"

void  
EqEMFieldXtal::SetChargeMomentumMass(G4double particleCharge, // e+ units
		                            G4double,
                                            G4double particleMass)
{
   fElectroMagCof =  eplus*particleCharge*c_light;
   //fElectroMagCof = electron_charge;
   fMassCof = particleMass*particleMass ; 
   /*
   G4cout<<"\nEqEMFIeldXtal::SetChargeMomentumMass: using normalToValley xform: \n\t"<<normalToValley.NetRotation().colX() <<"\n\t"<<normalToValley.NetRotation().colY()<<"\n\t"<<normalToValley.NetRotation().colZ();

   G4cout<<"\nEqEMFIeldXtal::SetChargeMomentumMass: using valleyToNormal xform: \n\t"<<valleyToNormal.NetRotation().colX() <<"\n\t"<<valleyToNormal.NetRotation().colY()<<"\n\t"<<valleyToNormal.NetRotation().colZ();

   G4cout<<"EqEMField::SetChargeMomentumMass: Using the affine xform, a "<<G4ThreeVector(1*mm,0,0)/mm<< " vector becomes: "<<valleyToNormal.TransformPoint(normalToValley.TransformPoint(G4ThreeVector(1*mm,0,0)))/mm;
   */
}



void
EqEMFieldXtal::EvaluateRhsGivenB(const G4double y[],
			                const G4double field[],
				              G4double dydx[] ) const
{ 
    G4ThreeVector pc = G4ThreeVector(y[3], y[4], y[5]);
    G4double pSquared = pc.mag2();
    G4double pModuleInverse  = 1.0/pc.mag();

    G4double Energy   = std::sqrt( pSquared + fMassCof );
    G4double cof1     = fElectroMagCof*pModuleInverse ;
    G4double cof2     = Energy/c_light ;

    dydx[0] = y[3]*pModuleInverse ;                         
    dydx[1] = y[4]*pModuleInverse ;                         
    dydx[2] = y[5]*pModuleInverse ;                        

    G4ThreeVector Efield = G4ThreeVector(field[3], field[4], field[5]);
    G4ThreeVector retForce = cof1*cof2*Efield;
    if (fElectroMagCof < 0)                                 //If electron
    {
      G4ThreeVector Efield_dir = Efield.unit();
      normalToValley.ApplyPointTransform(Efield_dir);
      G4ThreeVector rot_effect = G4ThreeVector(Efield_dir[0], 
				  Efield_dir[1], Efield_dir[2]);  
      for (int i = 0; i < 3; ++i) retForce[i] *= rot_effect[i];
      normalToValley.ApplyPointTransform(retForce); 
    }
    
    dydx[3] = retForce[0];
    dydx[4] = retForce[1]; 
    dydx[5] = retForce[2];  

    dydx[6] = 0.;//not used

    // Lab Time of flight
    G4double inverse_velocity = Energy * pModuleInverse / c_light;
    dydx[7] = inverse_velocity;
    return ;
    /*
   G4double pSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;

   G4double Energy   = std::sqrt( pSquared + fMassCof );
   G4double cof2     = Energy/c_light ;

   G4double pModuleInverse  = 1.0/std::sqrt(pSquared) ;

   G4double inverse_velocity = Energy * pModuleInverse / c_light;

   G4double cof1     = fElectroMagCof*pModuleInverse ;

   dydx[0] = y[3]*pModuleInverse ;                         
   dydx[1] = y[4]*pModuleInverse ;                         
   dydx[2] = y[5]*pModuleInverse ;                        

   G4ThreeVector retForce;
   retForce[0] = cof1*(cof2*field[3] + (y[4]*field[2] - y[5]*field[1]));
   retForce[1] = cof1*(cof2*field[4] + (y[5]*field[0] - y[3]*field[2]));
   retForce[2] = cof1*(cof2*field[5] + (y[3]*field[1] - y[4]*field[0]));
   
   dydx[3] = retForce[0];
   
   dydx[4] = retForce[1]; 
 
   dydx[5] = retForce[2];  

   dydx[6] = 0.;//not used

   // Lab Time of flight
   dydx[7] = inverse_velocity;
   //G4cout <<  "kmag: " <<  sqrt(pSquared/hbar_Planck)*m <<  G4endl;
   return ;
   */
/*
  double Y[6];
  Y[0]=posn.getX();
  Y[1]=posn.getY();
  Y[2]=posn.getZ();
  Y[3]=vel.getX();
  Y[4]=vel.getY();
  Y[5]=vel.getZ();


  double Field[6];  
  Field[0]=bField.getX();
  Field[1]=bField.getY();
  Field[2]=bField.getZ();
  Field[3]=eField.getX();
  Field[4]=eField.getY();
  Field[5]=eField.getZ();


  G4double pSquared = Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5] ;

  G4double Energy   = std::sqrt( pSquared + fMassCof );
  G4double cof2     = Energy/c_light ;


  G4double pModuleInverse  = 1.0/std::sqrt(pSquared) ;

  //  G4double inverse_velocity = Energy * c_light * pModuleInverse;
  G4double inverse_velocity = Energy * pModuleInverse / c_light;

  G4double cof1     = fElectroMagCof*pModuleInverse ;

  //  G4double vDotE = y[3]*Field[3] + y[4]*Field[4] + y[5]*Field[5] ;

  G4ThreeVector retVel = G4ThreeVector(Y[3]*pModuleInverse, Y[4]*pModuleInverse, Y[5]*pModuleInverse);
  G4ThreeVector retForce = G4ThreeVector(1,0,0);
  retForce.setX( 0.1*cof1*(cof2*Field[3] + (Y[4]*Field[2] - Y[5]*Field[1])) );
  retForce.setY( 0.1*cof1*(cof2*Field[4] + (Y[5]*Field[0] - Y[3]*Field[2])) );
  retForce.setZ( cof1*(cof2*Field[5] + (Y[3]*Field[1] - Y[4]*Field[0])) );

  //G4ThreeVector retVel = G4ThreeVector(Y[3],Y[4],Y[5])*pModuleInverse;
  //G4ThreeVector retForce = G4ThreeVector(cof2*Field[3] + (Y[4]*Field[2] - 
//Y[5]*Field[1]),
   //                                     cof2*Field[4] + (Y[5]*Field[0] - 
//Y[3]*Field[2]),
     //                                   cof2*Field[5] + (Y[3]*Field[1] - 
//Y[4]*Field[0]))
       //                                 *cof1;
                                        
  retVel=valleyToNormal.TransformPoint(retVel);
  retForce=valleyToNormal.TransformPoint(retForce);

  dydx[0] = retVel.getX() ;                         
  dydx[1] = retVel.getY() ;                         
  dydx[2] = retVel.getZ() ;                        

  dydx[3] = retForce.getX();  
  dydx[4] = retForce.getY();
  dydx[5] = retForce.getZ();
  
  
  dydx[6] = 0.;//not used

  // Lab Time of flight
  dydx[7] = inverse_velocity;
  return ;
  */
}

void EqEMFieldXtal::SetValleyTransform(G4AffineTransform xform){
  normalToValley = xform;
  valleyToNormal = xform.Inverse();
}
