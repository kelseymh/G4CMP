#include "globals.hh"

namespace G4CMP {
  template <class T> G4int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }
};
