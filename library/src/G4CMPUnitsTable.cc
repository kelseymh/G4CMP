/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
/// \file library/include/G4CMPUnitsTable.hh
/// \brief Define additional units with symbols for crystal applications
//

#include "G4CMPUnitsTable.hh"


// Create dummy class to trigger constructor below

G4ThreadLocal G4CMPUnitsTable* G4CMPUnitsTable::g4cmpUnitsTable = 0;

void G4CMPUnitsTable::Init() {
  if (!g4cmpUnitsTable) g4cmpUnitsTable = new G4CMPUnitsTable;
}

// Define special time units for scattering rate coefficients

G4CMPUnitsTable::G4CMPUnitsTable() {
  // Phonon scattering coeffient B [s^3]
  new G4UnitDefinition(     "second3","s3"   ,"Time cubed",second3);
  new G4UnitDefinition("millisecond3","ms3"  ,"Time cubed",millisecond3);
  new G4UnitDefinition("microsecond3","mus3" ,"Time cubed",microsecond3);

  // Phonon anharmonic decay coefficient A [s^4]
  new G4UnitDefinition(     "second4","s4"   ,"Time fourth",second4);
  new G4UnitDefinition("millisecond4","ms4"  ,"Time fourth",millisecond4);
  new G4UnitDefinition("microsecond4","mus4" ,"Time fourth",microsecond4);

  // Stiffness (pressure) and frequency units suitable for solid state physics
  new G4UnitDefinition("gigapascal", "GPa", "Pressure",  gigapascal);
  new G4UnitDefinition( "terahertz", "THz", "Frequency", terahertz);

  // Velocity (?!? Why aren't these already defined in Geant4 ?!?)
  new G4UnitDefinition(     "meters/second",  "m/s", "Velocity", m/s);
  new G4UnitDefinition( "kilometers/second", "km/s", "Velocity", km/s);
  new G4UnitDefinition("millimeters/second", "mm/s", "Velocity", mm/s);
  new G4UnitDefinition("centimeters/second", "cm/s", "Velocity", cm/s);
}


