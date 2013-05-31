
#include "TPhononSlow.hh"
#include "G4ParticleTable.hh"

TPhononSlow* TPhononSlow::theInstance = 0;
LatticeManager2* TPhononSlow::LM = NULL;

TPhononSlow*  TPhononSlow::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "TPhononSlow";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  // create particle
  //      
  //    Arguments for constructor are as follows 
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 
  //             shortlived      subType    anti_encoding
   anInstance = new G4ParticleDefinition(
                 name,         0.0*MeV,       0.0*MeV,         0.0,
                    0,               0,             0,
                    0,               0,             0,
             "phonon",               0,             0,         0,
                 true,             0.0,          NULL,
	        false,        "TPhononSlow",               0
	     );
  }
  theInstance = reinterpret_cast<TPhononSlow*>(anInstance);
  return theInstance;
}


TPhononSlow*  TPhononSlow::PhononDefinition() 
{
  return Definition();
}


void TPhononSlow::setLatticeManager(LatticeManager2* Lat){
  LM=Lat;
}

LatticeManager2* TPhononSlow::getLatticeManager(){
  return LM;
}
