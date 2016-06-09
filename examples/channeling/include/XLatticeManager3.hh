/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/XLatticeManager3.hh
/// \brief Definition of the XLatticeManager3 class
//
// $Id: b40ddc55e39064f71b49611db914c9b55adcceea $
//
#ifndef XLatticeManager3_h
#define XLatticeManager3_h 1


#define MAXLAT 10 //maximum number of different crystal lattices supported

#include "G4ThreeVector.hh"
#include "XPhysicalLattice.hh"

class XLatticeManager3
{
private:
  static XLatticeManager3* LM;

protected:
  XLatticeManager3();
  ~XLatticeManager3();

  XPhysicalLattice* fLatticeList[MAXLAT]; 
  int fTotalLattices;

public:

  static XLatticeManager3* GetXLatticeManager(); 

  XPhysicalLattice* GetXPhysicalLattice(G4VPhysicalVolume*);
  bool RegisterLattice(XPhysicalLattice*);
  bool HasLattice(G4VPhysicalVolume*);
  double MapKtoV(G4VPhysicalVolume*,int,const G4ThreeVector &);
  G4ThreeVector MapKtoVDir(G4VPhysicalVolume*,int,const G4ThreeVector &);

};



#endif
