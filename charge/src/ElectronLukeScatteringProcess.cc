#include "ElectronLukeScatteringProcess.hh"
#include "DriftingElectron.hh"
#include "DriftingHole.hh"

#include "LPhonon.hh"
#include "PhononTrackInformation.hh"
#include "LatticeManager2.hh"
#include "PhysicalLattice.hh"

#include "G4FieldManager.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include <fstream>
#include <iostream>
#include "G4TransportationManager.hh"
#include "DriftingElectronTrackInformation.hh"
#include "G4Field.hh"
#include "Tst1EMField.hh"

ElectronLukeScatteringProcess::ElectronLukeScatteringProcess(G4VProcess* sLim)
 : G4VDiscreteProcess("ElectronLukeScattering"), velLong(5324.2077*m/s), 
		    l0(257e-6*m), me(electron_mass_c2/c_squared),
		    mc(.1185*electron_mass_c2/c_squared), stepLimiter(sLim)
{
    ksound = mc*velLong/hbar_Planck;
    G4cout << "ElectronLukeScatteringProcess::Constructor:ksound_Hole = "
	   << ksound*m << " /m" << G4endl;
    
    if(verboseLevel > 1)
    {
	G4cout << GetProcessName() << " is created." << G4endl;
    }
}

ElectronLukeScatteringProcess::~ElectronLukeScatteringProcess()
{ ; }

ElectronLukeScatteringProcess::ElectronLukeScatteringProcess(ElectronLukeScatteringProcess& right)
  : G4VDiscreteProcess(right)
{ ; }

G4double ElectronLukeScatteringProcess::GetMeanFreePath(const G4Track& aTrack, 
					G4double, G4ForceCondition* condition)
{
    G4RotationMatrix trix;
    int valley = 
	((DriftingElectronTrackInformation*) aTrack.GetUserInformation())
								->getValley();
    switch(valley) {
    case 1:
        trix = G4RotationMatrix(
		    G4Rep3x3( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
                             -1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
                             -1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
        break;
    case 2:
        trix = G4RotationMatrix(
		    G4Rep3x3(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0),
                             -1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
                              1.0/sqrt(6.0), -1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
        break;
    case 3:
        trix = G4RotationMatrix(
		    G4Rep3x3(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
                              1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0,
                              1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
        break;
    case 4:
        trix = G4RotationMatrix(
		    G4Rep3x3( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),
                              1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0,
                             -1.0/sqrt(6.0),  1.0/sqrt(6.0),  sqrt(2.0/3.0) ));
        break;
    }

    ValleyToNormal= G4AffineTransform(trix);
    NormalToValley= G4AffineTransform(trix).Inverse();
    T = G4ThreeVector(sqrt(.1185/1.588), sqrt(.1185/.081), sqrt(.1185/.081));
    
//    G4FieldManager* fMan =
//        G4TransportationManager::GetTransportationManager()->GetFieldManager();
//    const G4Field* field = fMan->GetDetectorField();

//    G4ThreeVector posVec = aTrack.GetPosition();
//    G4double position[4] = {posVec[0],posVec[1],posVec[2],0};
//    G4double fieldVal[6];

//    field->GetFieldValue(position,fieldVal);
//    G4ThreeVector Efield = G4ThreeVector(fieldVal[3], fieldVal[4], fieldVal[5]);
//    G4ThreeVector Efield_valley = NormalToValley.TransformPoint(Efield);
//    G4ThreeVector Efield_HV = G4ThreeVector( Efield_valley[0]*T[0],
//               Efield_valley[1]*T[1],
//               Efield_valley[2]*T[2]);

/*    mInv = trix.inverse()*G4Rep3x3(1/1.588/me,   0.0    , 0.0,
                                   0.0       , 1/.081/me, 0.0,
                                   0.0       ,   0.0    , 1/.081/me)*trix*/;

    G4ThreeVector k = aTrack.GetMomentum()/hbarc;
    G4ThreeVector v = mInv*k*hbar_Planck;
    G4ThreeVector k_valley = NormalToValley.TransformPoint(k);
    G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0],
                                        k_valley[1]*T[1],
                                        k_valley[2]*T[2]);

    *condition = Forced;

    G4double kmag = k_HV.mag();

    if(kmag<=ksound)
    {
        return DBL_MAX;
    }

    G4double tau0 =  1.0 / (
                         velLong / (3*l0)
                         * (kmag / ksound) * (kmag/ksound)
                         * ((1- ksound /kmag))
                         * ((1- ksound /kmag))
                         * ((1- ksound /kmag))
                     );

    G4double mfp0 = tau0*v.mag();
    return mfp0;
    /*
    G4cout << mfp0/m << G4endl;

    G4ThreeVector k1, k1_HV, x1, v1;
    G4ThreeVector v0 = v;
    G4ThreeVector x0 = posVec;
    G4double mfp1, tau1;
    G4double dmfp = 0;
    for (int i=0; i<100; i++)
    {
        x1 = x0 + mfp0*v0/v0.mag();
        k1_HV =  k_HV + electron_charge/hbar_Planck*Efield_HV*mfp0/v0.mag()/100;
        v1 = v0 + mInv*Efield*electron_charge*mfp0/v0.mag()/100;

        //G4cout << (k1_HV.mag() - k_HV.mag())*m << G4endl;

        tau1 =  1.0 / (
                    velLong / (3*l0)
                    * (k1_HV.mag()/ ksound)
                    * ((1- ksound /k1_HV.mag()))
                    * ((1- ksound /k1_HV.mag()))
                    * ((1- ksound /k1_HV.mag()))
                );

        //G4cout << "mfp0 = " << mfp0/m << G4endl;
        //G4cout << "mfp1 = " << tau1*v1.mag()/m << G4endl;
        //G4cout << "mfp2 = " << tau1*v1.mag()/m  + mfp0/10/m<< G4endl;
        dmfp += mfp0/100;
        mfp0 = mfp0/100 + tau1*v1.mag();
        x0 = x1;
        k_HV = k1_HV;
        v0 = v1;
    }
    //dmfp += tau1*v1.mag();
//   for (int i=0; i<100; i++)
//   {
//   //x1 = x0 + mfp0*v0/v0.mag();
//   k1_HV =  k_HV + electron_charge/hbar_Planck*Efield_HV*mfp0/v0.mag()/100;
//   v1 = v0 + mInv*Efield*electron_charge*mfp0/v0.mag()/100;
//
//   G4cout << (k1_HV.mag() - k_HV.mag())*m << G4endl;
//
//   tau1 =  1.0 / (
// 		       velLong / (3*l0Hole)
// 		       * (k1_HV.mag()/ ksound_Hole)
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       * ((1- ksound_Hole /k1_HV.mag()))
// 		       );
//
//     G4cout << "mfp0 = " << mfp0/m << G4endl;
//     G4cout << "mfp1 = " << tau1*v1.mag()/m << G4endl;
//     G4cout << "mfp2 = " << tau1*v1.mag()/m  + mfp0/10/m<< G4endl;
//     mfp0 = mfp0/100 + tau1*v1.mag();
//   }
    G4cout << mfp0/m << G4endl;
    return dmfp;*/


}

G4VParticleChange* ElectronLukeScatteringProcess::PostStepDoIt(
						    const G4Track& aTrack,
						    const G4Step& aStep )
{
    aParticleChange.Initialize(aTrack);
    G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();

    G4ThreeVector k = postStepPoint->GetMomentum()/hbarc;
    G4ThreeVector k_valley = NormalToValley.TransformPoint(k);
    G4ThreeVector k_HV =  G4ThreeVector(k_valley[0]*T[0],
                                        k_valley[1]*T[1],
                                        k_valley[2]*T[2]);

    G4double kmag = k_HV.mag();

    //Do nothing other than re-calculate mfp when step limit reached or leaving
    //volume
    if((postStepPoint->GetProcessDefinedStep()==stepLimiter)
            ||(postStepPoint->GetStepStatus()==fGeomBoundary)
            || (kmag <= ksound))
    {
        return &aParticleChange;
    }

    G4double theta_phonon = MakeTheta(kmag, ksound);
    G4double theta_charge = 
			acos( (kmag*kmag - 2*ksound*(kmag*cos(theta_phonon)
								    - ksound)
                               - 2 * (kmag*cos(theta_phonon) - ksound)
				   * (kmag*cos(theta_phonon) - ksound) )/
                              kmag/ (sqrt(kmag*kmag - 4*ksound
                                   * (kmag*cos(theta_phonon) - ksound) ) ) );

//   std::ofstream charge("theta_charge", std::ios::app);
//   std::ofstream phonon("theta_phonon", std::ios::app);
//   charge << theta_charge << G4endl;
//   phonon << theta_phonon << G4endl;
//   charge.close();
//   phonon.close();

    G4double q = 2*(kmag*cos(theta_phonon)-ksound);
    //k_HV.setMag(sqrt(k_HV.mag2() + q*q - 2*kmag*q*cos(theta_phonon)));
    k_HV.setMag(sqrt(kmag*kmag-2*mc*q*velLong/hbar_Planck));
    G4ThreeVector kdir = k_HV.unit();
    k_HV.rotate(kdir.orthogonal(), theta_charge);

    G4double phi_charge =  G4UniformRand()*2*pi;
    k_HV.rotate(kdir, phi_charge);

    G4ThreeVector p_new = hbar_Planck*k_HV;
    p_new[0] /= T[0];
    p_new[1] /= T[1];
    p_new[2] /= T[2];
    ValleyToNormal.ApplyPointTransform(p_new);

    G4ThreeVector phononq = q*k.unit().rotate(k.unit(), theta_phonon);
//   std::ofstream charge("theta_charge", std::ios::app);
    std::ofstream phonon("theta_phonon", std::ios::app);
//   charge << theta_charge << G4endl;
    phonon << phononq.getZ()/phononq.mag()<< G4endl;
//   charge.close();
    phonon.close();
    /*
      std::ofstream epositions;
      epositions.open("qenergies.txt", std::ofstream::app);
    files
      epositions << aTrack.GetKineticEnergy()/eV - Energy << "\n";
      epositions.close();
      */


//    std::ofstream energy;
//    energy.open("energy", std::ofstream::app);
//
//    energy<< aTrack.GetKineticEnergy()/eV << "\n";
//    energy.close();

//     std::ofstream qwave;
//     qwave.open("q", std::ofstream::app);
//     qwave << hbar_Planck*velLong*q/eV << "\n";
//     qwave.close();
//
//   std::ofstream kwave;
//   kwave.open("k", std::ofstream::app);
//
//   kwave << oldE/eV - Energy/eV << "\n";
//   kwave.close();
   std::ofstream espec("energy_spectrum_phonon", std::ios::app);
   //espec << (-hbar_Planck*hbar_Planck*k_HV.mag2()/2.0/mc + p_new.mag2()/2.0/mc)/h_Planck/hertz << G4endl;
   espec << (k_HV.mag2()*hbarc_squared/c_squared/2/mc)/eV*1000<< G4endl;
   //espec << (p_new.mag2()/2/mc)/eV*1000<< G4endl;
   espec.close();


    aParticleChange.ProposeMomentumDirection(p_new.unit());
    aParticleChange.ProposeEnergy(p_new.mag2()/2/mc);
    ResetNumberOfInteractionLengthLeft();
    return &aParticleChange;
}

G4double ElectronLukeScatteringProcess::MakeTheta(G4double& k, G4double& ks)
{
    //Analytical method to compute theta
    G4double u = G4UniformRand();
    G4double base = -(u-1)+3*(u-1)*(ks/k)
		    -3*(u-1)*(ks/k)*(ks/k)
		    +(u-1)*(ks/k)*(ks/k)*(ks/k);
    G4double exponent = 1.0/3.0;
    G4double operand = ks/k+pow(base, exponent);

    if(operand>1.0) {
        // G4cout<<"\nTruncating operand from"<<operand<<" to 1.0";
        operand=1.0;
    }
    if(operand<0.0) G4cout<<"\noperand out of range: operand = "<<operand;
	
    G4double theta = acos(operand);

    if(acos(ks/k)<theta) G4cout<<"\n  THETA OUT OF BOUNDS!!! (theta>acos(ks/k)) theta:"<<theta<<" acos(ks/k):"<<acos(ks/k);
    if(pi/2<=theta) G4cout<<"\n THETA OUT OF BOUNDS!!! (pi/2 < theta)";

//   //Rejection method for determining theta
//   G4double theta;
//   G4double pValue=2;
//   G4double pDensity = 1;
//
//   G4double thetaMax;
//   G4double pValMax;
//   thetaMax=acos(ks/k);
//   pValMax=(1-(ks/k))*(1-(ks/k));
//   bool first = false;
//
//   while(pValue>pDensity){
//     pValue=G4UniformRand();//*pValMax;//  *(1+2*ks/k+(ks/k)*(ks/k));
//     theta=G4UniformRand()*thetaMax;
//     pDensity = (cos(theta)-(ks/k))*(cos(theta)-(ks/k))*sin(theta);//need to multiply by unit 's' to make it dimensionless
//     if(pDensity>pValMax) G4cout<<"\nLukeScattering::PostStepDoIt: Error: pDensity should never exceed pValMax "<<pValMax;


    //    if(!first){
    //G4cout<<"\nLukeScatteringProcess::MakeTheta: pDensity calculated as:
//"<<pDensity;
    //   G4cout<<"\n\tpValue: "<<pValue;
    // G4cout<<"\n\ttheta: "<<theta/rad;
    // G4cout<<"\n\tk: "<<k*m;
    //  G4cout<<"\n\tks: "<<ks*m;
    //  G4cout<<endl;
    //}
    //first=false;

//  }
  return theta;
}

G4double ElectronLukeScatteringProcess::MakePhi(G4double& k,G4double& ks, 
						G4double& theta)
{
    G4double phi = acos((k*k - 2*ks*(k*cos(theta)-ks)-2*(k*cos(theta)-ks)*(k*cos(theta)-ks))/(k*sqrt(k*k-4*ks*(k*cos(theta)-ks))));

    return phi;
}

G4bool ElectronLukeScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
    return((&aPD==DriftingElectron::Definition()));
}
