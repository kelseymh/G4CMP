/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file HeterostructureActionInitialization.hh
/// \brief Definition of the HeterostructureActionInitialization class

// 20260815  Selby Q. Dang (Stanford/SLAC)

#ifndef HeterostructureActionInitialization_hh
#define HeterostructureActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class HeterostructureActionInitialization : public G4VUserActionInitialization {
public:
  HeterostructureActionInitialization() {;}
  virtual ~HeterostructureActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* HeterostructureActionInitialization_hh */
