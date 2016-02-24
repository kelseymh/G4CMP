/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4CMPDriftHole.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4CMPDriftHole* G4CMPDriftHole::theInstance = 0;

G4CMPDriftHole* G4CMPDriftHole::Definition(){  

  if (theInstance !=0) return theInstance;

  // search in particle table
  const G4String name = "G4CMPDriftHole";
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
		    name,         0.35*electron_mass_c2,       0.0*MeV,         +1.*eplus,
                    1,               0,             0,
                    0,               0,             0,
             "lepton",               0,             0,         0,
                 true,             0.0,          NULL,
	        false,        "G4CMPDriftHole",   0
		 );
  }

  theInstance = reinterpret_cast<G4CMPDriftHole*>(anInstance);
  return theInstance;


}

G4CMPDriftHole* G4CMPDriftHole::G4CMPDriftHoleDefinition(){
  return Definition();
}
