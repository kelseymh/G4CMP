#include "QPDownconversion.hh"

#include "TPhononFast.hh"
#include "TPhononSlow.hh"
#include "LPhonon.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "PhononTrackInformation.hh"

#include <cmath>
#include <stack>



QPDownconversion::QPDownconversion(){

writer.open("timing2.ssv",fstream::in | fstream::out | fstream::ate);

}

QPDownconversion::~QPDownconversion(){

  writer.close();

}

void QPDownconversion::downconvert(const G4Track& aTrack, G4VParticleChange* aParticleChange, G4ThreeVector direction){

  //int nos = 100; //number of secondaries
  //aParticleChange->SetNumberOfSecondaries(nos);
  G4double alThick = 300*nm;
  G4double alInteract = 720*nm;
  G4double alGap = 2*170e-6*eV;
  G4double edep = 0*eV;

  //  G4cout<<"\nQPDownconversion::downconvert: ...has been called with particle energy "<<aTrack.GetTotalEnergy()/eV<<" eV.";
  //  G4cout<<"\nQPDownconversion::downconvert: Likeliehood of phonon-QP interaction: "<<exp(-2*2*alThick/alInteract);

  //Does a phonon/QP interaction occur?
  bool interact = (G4UniformRand()<exp(-2*2*alThick/alInteract));
  if(interact){
    //    G4cout<<"\nQPDownconverison::downconvert: phonon interaction with particle ocurring...";
    stack <G4DynamicParticle*> phonons;
    stack <G4double> QPs;//only need doubles since only attribute of interest is energy

    phonons.push(new G4DynamicParticle(*(aTrack.GetDynamicParticle())));
    
    while((phonons.size()>0)|(QPs.size()>0)){
      //If there are phonons, break cooper paris
      if(phonons.size()>0){
	G4double pE = phonons.top()->GetTotalEnergy();
	if(pE>2*alGap){
	  QPs.push(QPEnergy(pE));
	  QPs.push(pE-QPs.top());
	}
	delete(phonons.top());
	phonons.pop();
      }
      //if there are quasi particles, emit phonons
      if(QPs.size()>0){
	if(QPs.top()>3*alGap){
	  phonons.push(new G4DynamicParticle(LPhonon::Definition(), direction, QPs.top()-PhEnergy(QPs.top())));
	  QPs.top()=QPs.top()-phonons.top()->GetKineticEnergy();
	}
	else{
	  //if QP energy is low, delete QP and increment edep by algap
	  edep=edep+QPs.top();
	  writer<<aTrack.GetGlobalTime()<<"\n";
	  QPs.pop();
	}
      }
      //if there are phonons, see if they escape the Al
      //currently all phonons escape
      while(phonons.size()>0){
	if(G4UniformRand()<0.4728){
	  //	aParticleChange->AddSecondary(new G4DynamicParticle());
	  G4Track* sec = new G4Track(phonons.top(), aTrack.GetGlobalTime(), aTrack.GetPosition());
	  sec->SetUserInformation(new PhononTrackInformation(direction));
	  //nos++;
	  //aParticleChange->SetNumberOfSecondaries(nos);
	  aParticleChange->AddSecondary(sec);
	  //G4cout<<"\nQPDownconversion::downconvert: emitting phonon back into crystal with energy "<<phonons.top()->GetTotalEnergy()/alGap;
	  //	G4cout<<"\nQPDownconversion::downconvert: Adding secondary...\n";
	  phonons.pop();
	}
	
      }
    }
    //if there was an interaction, record energy deposition and delete original track
    aParticleChange->ProposeNonIonizingEnergyDeposit(edep);
    aParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  //if there was no interaction, just return and leave reflected phonon unchanged  
}

G4double QPDownconversion::PhEnergy(G4double QPEi){
  G4double gap = 2*170e-6*eV;
  G4double QPE = 0.5*QPEi;
  G4double P=1000;
  G4double llimit = gap;
  G4double ulimit = QPEi;
  while((QPEi-QPE)*(QPEi-QPE)*(1-(gap*gap)/(QPEi*(QPEi-(QPEi-QPE))))*((QPEi-(QPEi-QPE))/sqrt((QPEi-(QPEi-QPE))*(QPEi-(QPEi-QPE))-gap*gap))<P){
    QPE=llimit+(ulimit-llimit)*G4UniformRand();
    P=G4UniformRand()*100;
  }
  return QPE;
}

G4double QPDownconversion::QPEnergy(G4double phononE){
  G4double gap = 2*170e-6*eV;
  G4double QPE = 0.5*phononE;
  G4double P=1000;
  G4double llimit = gap;
  G4double ulimit = QPE-gap;
  
  while( (1/sqrt(QPE*QPE-gap*gap))*(QPE*(phononE-QPE)+gap*gap)/sqrt((phononE-QPE)*(phononE-QPE)-gap*gap)<P){
    P=G4UniformRand()*100;
    QPE=llimit+(ulimit-llimit)*G4UniformRand();
  }
  return QPE;
}
