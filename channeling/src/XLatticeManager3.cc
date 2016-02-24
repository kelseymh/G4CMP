/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/XLatticeManager3.cc
/// \brief Implementation of the XLatticeManager3 class
//
// $Id: 48267b45e4b280b360dbd946c7eeeae4d4462677 $
//

#include "XLatticeManager3.hh"
#include "G4VPhysicalVolume.hh"

//int XLatticeManager3::fTotalLattices = 0;
XLatticeManager3* XLatticeManager3::LM;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XLatticeManager3::XLatticeManager3()
{
  fTotalLattices = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLatticeManager3::~XLatticeManager3()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLatticeManager3* XLatticeManager3::GetXLatticeManager(){

  //if no lattice manager exists, create one.
  if(!LM) LM = new XLatticeManager3();

  //return pointer to single existing lattice manager
  return LM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool XLatticeManager3::RegisterLattice(XPhysicalLattice* Lat){
  if(fTotalLattices<MAXLAT){
    fTotalLattices++;
    //using "fTotalLattices-1" so that first lattice corresponds to index 0
    fLatticeList[fTotalLattices-1]=Lat;  
    G4cout<<"\nXLatticeManager3::registerLattice: Registering Lattice. Total number of lattices:"<<fTotalLattices<<"\n"<<endl;

    return true; 
  }
  
  G4cout<<"\nXLatticeManager::RegisterLattice(XPhysicalLattice*): Maximum number of lattices MAXLAT exceeded."<<endl;
 
  return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* XLatticeManager3::GetXPhysicalLattice(G4VPhysicalVolume* Vol){
//returns a pointer to the PhysicalLattice associated with Vol

    for(int counter=0;counter<fTotalLattices;counter++){
      if(fLatticeList[counter]->GetVolume()==Vol) {
        return fLatticeList[counter]; //found matching lattice
      }      
    }
    return fLatticeList[0];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool XLatticeManager3::HasLattice(G4VPhysicalVolume* Vol){
  //return true if Vol has a physical lattice

    for(int counter=0;counter<fTotalLattices;counter++){
      if(fLatticeList[counter]->GetVolume()==Vol) {
        return true; //found matching lattice
      }      
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

double XLatticeManager3::MapKtoV(G4VPhysicalVolume* Vol, int polarizationState, const G4ThreeVector & k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon velocity in m/s

  if((Vol==NULL)&&(fTotalLattices>0)) return fLatticeList[0]->MapKtoV(polarizationState, k);
  for(int counter=0;counter<fTotalLattices;counter++){
    if(fLatticeList[counter]->GetVolume()==Vol) {
      return fLatticeList[counter]->MapKtoV(polarizationState, k);
    }
  }
  G4cout<<"\nXLatticeManager::MapKtoV: Found no matching lattices for " <<Vol->GetName()<<". Total number of lattices is "<<fTotalLattices<<endl;
return 300;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XLatticeManager3::MapKtoVDir(G4VPhysicalVolume* Vol, int polarizationState, const G4ThreeVector & k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon propagation direction as dimensionless unit vector

  if((Vol==NULL)&&(fTotalLattices>0)) return fLatticeList[0]->MapKtoVDir(polarizationState, k);
  for(int counter=0;counter<fTotalLattices;counter++){
    if(fLatticeList[counter]->GetVolume()==Vol) {
      if(counter!=0) G4cout<<"\nLattiveManager2::MapKtoV:returning group velocity from lattise position: "<<counter;
      return fLatticeList[counter]->MapKtoVDir(polarizationState, k);
    }
  }
  //G4cout<<"\nXLatticeManager::MapKtoVDir: Found no matching lattices for " <<Vol->GetName()<<". Total number of lattices is "<<fTotalLattices<<endl;
return G4ThreeVector(1,0,0);
}

