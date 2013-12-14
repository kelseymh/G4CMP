
#include "LPhonon.hh"
#include "G4ParticleTable.hh"

LPhonon* LPhonon::theInstance = 0;

LPhonon*  LPhonon::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "LPhonon";
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
	        false,        "LPhonon",               0
	     );
  }
  theInstance = reinterpret_cast<LPhonon*>(anInstance);
  return theInstance;
}

LPhonon*  LPhonon::PhononDefinition() 
{
  return Definition();
}


