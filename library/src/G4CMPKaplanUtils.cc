/***********************************************************************\
 *  * This software is licensed under the terms of the GNU General Public *
 *   * License version 3 or later. See G4CMP/LICENSE for the full license. *
 *   \***********************************************************************/

/// \file library/src/G4CMPKaplanUtils.cc
//
//


#include "globals.hh"
#include "G4CMPKaplanUtils.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <numeric>

G4double G4CMPKaplanUtils::QPEnergyRand(G4double Energy) const {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  // Instead, we'll do a rejection method.
  // 
  // PDF(E') = (E'*(Energy - E') + gapEnergy*gapEnergy)
  //          /
  //           sqrt((E'*E' - gapEnergy*gapEnergy) *
  //              ((Energy - E')*(Energy - E') - gapEnergy*gapEnergy));
  // The shape of the PDF is like a U, so the max values are at the endpoints:
  // E' = gapEnergy and E' = Energy - gapEnergy
  
  // Add buffer so first/last bins don't give zero denominator in pdfSum


  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + (Energy-2.*gapEnergy)/BUFF;
  G4double xmax = gapEnergy + (Energy-2.*gapEnergy)*(BUFF-1.)/BUFF;
  G4double ymax = QPEnergyPDF(Energy, xmin);

  G4double xtest=0., ytest=ymax;
  do {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  } while (ytest > QPEnergyPDF(Energy, xtest));

  return xtest;
}

G4double G4CMPKaplanUtils::QPEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  return ( (x*(E-x) + gapsq) / sqrt((x*x-gapsq) * ((E-x)*(E-x)-gapsq)) );
}

G4double G4CMPKaplanUtils::PhononEnergyRand(G4double Energy) const {
  // PDF is not integrable, so we can't do an inverse transform sampling.
  //   Instead, we'll do a rejection method.
  //  /
  //   PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
  //             /
  //             sqrt((E'*E' - gapEnergy*gapEnergy);
  //

  const G4double BUFF = 1000.;
  G4double xmin = gapEnergy + gapEnergy/BUFF;
  G4double xmax = Energy;
  G4double ymax = PhononEnergyPDF(Energy, xmin);

  G4double xtest=0., ytest=ymax;
  do {
    ytest = G4UniformRand()*ymax;
    xtest = G4UniformRand()*(xmax-xmin) + xmin;
  } while (ytest > PhononEnergyPDF(Energy, xtest));

  return Energy-xtest;
}

G4double G4CMPKaplanUtils::PhononEnergyPDF(G4double E, G4double x) const {
  const G4double gapsq = gapEnergy*gapEnergy;
  return ( x*(E-x)*(E-x) * (x-gapsq/E) / sqrt(x*x - gapsq) );
}
