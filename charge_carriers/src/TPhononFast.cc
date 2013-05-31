
#include "TPhononFast.hh"
#include "G4ParticleTable.hh"

TPhononFast* TPhononFast::theInstance = 0;
LatticeManager2* TPhononFast::LM = NULL;

TPhononFast*  TPhononFast::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "TPhononFast";
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
	        false,        "TPhononFast",               0
	     );
  }
  theInstance = reinterpret_cast<TPhononFast*>(anInstance);
  return theInstance;
}


TPhononFast*  TPhononFast::PhononDefinition() 
{
  return Definition();
}


void TPhononFast::setLatticeManager(LatticeManager2* Lat){
  LM=Lat;
}

LatticeManager2* TPhononFast::getLatticeManager(){
  return LM;
}
