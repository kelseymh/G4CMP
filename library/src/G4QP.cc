/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file particles/src/G4QP.cc
/// \brief Implementation of the G4QP class
// This class deals with quasiparticles in the superconducting layer

#include "G4QP.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

G4QP* G4QP::theInstance = 0;


G4QP* G4QP::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "QP";
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
                 name,         100*eV,       0.0*MeV,         0.0,
                    1,               0,             0,
                    0,               0,             0,
             "QP",         0,             0,         0,
                 true,             0.0,          NULL,
                false,        "QP",           0
             );
  }
  theInstance = reinterpret_cast<G4QP*>(anInstance);
  return theInstance;
}


G4QP* G4QP::QPDefinition() {
  return Definition();
}


