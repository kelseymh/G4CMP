/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/phonons/src/G4PhononTransFast.cc
/// \brief Implementation of the G4PhononTransFast class
//
// $Id$
//

#include "G4PhononTransFast.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

G4PhononTransFast* G4PhononTransFast::theInstance = 0;

G4PhononTransFast*  G4PhononTransFast::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "phononTF";
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
                false,        "phononTF",       0
             );
  }
  theInstance = reinterpret_cast<G4PhononTransFast*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



G4PhononTransFast*  G4PhononTransFast::PhononDefinition() 
{
  return Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




