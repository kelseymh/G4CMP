/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/*
 *  \file electromagnetic/TestEm7/include/G4LindhardPartition.hh
 *  \brief Definition of the G4LindhardPartition class
 *
 *  Created by Marcus Mendenhall on 1/14/08.
 *  2008 Vanderbilt University, Nashville, TN, USA.
 *
 */

// $Id: 052e93dece41ae493dbef27183b4857142498437 $
//

#include "globals.hh"

class G4Material;

class G4VNIELPartition 
{
public:
        G4VNIELPartition() { }
        virtual ~G4VNIELPartition() { }
        
        // return the fraction of the specified energy which will be deposited as NIEL
        // if an incoming particle with z1, a1 is stopped in the specified material
        // a1 is in atomic mass units, energy in native G4 energy units.
        virtual G4double PartitionNIEL(
                G4int z1, G4double a1, const G4Material *material, G4double energy
        ) const =0;
};

class G4LindhardRobinsonPartition : public G4VNIELPartition
{
public:
        G4LindhardRobinsonPartition();
        virtual ~G4LindhardRobinsonPartition() { }
        
        virtual G4double PartitionNIEL(
                G4int z1, G4double a1, const G4Material *material, G4double energy
        ) const ;
        
        G4double z23[120];
        size_t   max_z;
};

