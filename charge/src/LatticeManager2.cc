
#include "LatticeManager2.hh"
#include "G4VPhysicalVolume.hh"

int LatticeManager2::totalLattices = 0;
PhysicalLattice* LatticeManager2::LatticeList[MAXLAT];

LatticeManager2::LatticeManager2()
{;}

LatticeManager2::~LatticeManager2()
{;}

bool LatticeManager2::registerLattice(PhysicalLattice* Lat){
  if(totalLattices<MAXLAT){
    totalLattices++;
    LatticeList[totalLattices-1]=Lat;  //using "totalLattices-1" so that first lattice corresponds to index 0
    G4cout<<"\nLatticeManager2::registerLattice: Registering Lattice. Total number of lattices:"<<totalLattices<<"\n"<<endl;
    return true; 
  }
  
  G4cout<<"\nLatticeManager::RegisterLattice(PhysicalLattice*): Maximum number of lattices MAXLAT exceeded."<<endl;
 
  return false;
}

PhysicalLattice* LatticeManager2::getPhysicalLattice(G4VPhysicalVolume* Vol){
    for(int counter=0;counter<totalLattices;counter++){
      if(LatticeList[counter]->getVolume()==Vol) {
	return LatticeList[counter]; //found matching lattice
      }      
    }
    return LatticeList[0];
}

bool LatticeManager2::hasLattice(G4VPhysicalVolume* Vol){
    for(int counter=0;counter<totalLattices;counter++){
      if(LatticeList[counter]->getVolume()==Vol) {
	return true; //found matching lattice
      }      
    }
    return false;
}

double LatticeManager2::mapKtoV(G4VPhysicalVolume* Vol, int polarizationState, const G4ThreeVector & k)
{

  if((Vol==NULL)&&(totalLattices>0)) return LatticeList[0]->mapKtoV(polarizationState, k);
  for(int counter=0;counter<totalLattices;counter++){
    if(LatticeList[counter]->getVolume()==Vol) {
      return LatticeList[counter]->mapKtoV(polarizationState, k);
    }
  }
  G4cout<<"\nLatticeManager::mapKtoV: Found no matching lattices for " <<Vol->GetName()<<". Total number of lattices is "<<totalLattices<<endl;
return 300;
}

G4ThreeVector LatticeManager2::mapKtoVDir(G4VPhysicalVolume* Vol, int polarizationState, const G4ThreeVector & k)
{

  if((Vol==NULL)&&(totalLattices>0)) return LatticeList[0]->mapKtoVDir(polarizationState, k);
  for(int counter=0;counter<totalLattices;counter++){
    if(LatticeList[counter]->getVolume()==Vol) {
      if(counter!=0) G4cout<<"\nLattiveManager2::MapKtoV:returning group velocity from lattise position: "<<counter;
      return LatticeList[counter]->mapKtoVDir(polarizationState, k);
    }
  }
  //G4cout<<"\nLatticeManager::mapKtoVDir: Found no matching lattices for " <<Vol->GetName()<<". Total number of lattices is "<<totalLattices<<endl;
return G4ThreeVector(1,0,0);
}

