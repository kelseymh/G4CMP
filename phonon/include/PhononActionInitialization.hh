// $Id: e58a61fedbb99b167e16dafebc9c8664ae0c7b94 $

#ifndef PhononActionInitialization_hh
#define PhononActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class PhononActionInitialization : public G4VUserActionInitialization {
public:
  PhononActionInitialization() {;}
  virtual ~PhononActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* PhononActionInitialization_hh */
