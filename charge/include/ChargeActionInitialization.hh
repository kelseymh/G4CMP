// $Id: 10d985fa839a054ee76ee8f3176a55e43ee840ce $
//

#ifndef ChargeActionInitialization_hh
#define ChargeActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class ChargeActionInitialization : public G4VUserActionInitialization {
public:
  ChargeActionInitialization() {;}
  virtual ~ChargeActionInitialization() {;}
  virtual void Build() const;
};

#endif	/* ChargeActionInitialization_hh */
