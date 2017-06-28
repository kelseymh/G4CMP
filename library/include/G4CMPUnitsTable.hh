#ifndef G4CMPUnitsTable_hh
#define G4CMPUnitsTable_hh 1
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
/// \file library/include/G4CMPUnitsTable.hh
/// \brief Define additional units with symbols for crystal applications
//
// 20160802 Use hep_pascal for pressure (Windows compatibility)
// 20170525 Drop unnecessary empty destructor

#include "G4Types.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


// Empty class to trigger defining units at initialization time

class G4CMPUnitsTable {
public:
  static void Init();		// Provides thread-by-thread initialization

private:
  static G4ThreadLocal G4CMPUnitsTable* g4cmpUnitsTable;
  G4CMPUnitsTable();
};

// Define global variables for additional units, with automatic "using"

static const double second3      = second*second*second;
static const double millisecond3 = millisecond*millisecond*millisecond;
static const double microsecond3 = microsecond*microsecond*microsecond;

static const double second4      = second3*second;
static const double millisecond4 = millisecond3*millisecond;
static const double microsecond4 = microsecond3*microsecond;

static const double gigapascal = 1e9*hep_pascal;
static const double terahertz  = 1e12*hertz;

static const double s3   = second3;
static const double ms3  = millisecond3;
static const double mus3 = microsecond3;

static const double s4   = second4;
static const double ms4  = millisecond4;
static const double mus4 = microsecond4;

static const double GPa  = gigapascal;
static const double THz  = terahertz;

#endif	/* G4CMPUnitsTable_hh */
