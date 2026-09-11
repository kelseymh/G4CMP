#ifndef G4CMPInterfaceTableGenerator_hh
#define G4CMPInterfaceTableGenerator_hh 1

#include "G4CMPAnisotropicInterfaceSolver.hh"
#include "G4CMPInterfaceTable.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include <cstddef>
#include <string>

#include "globals.hh"

class G4LatticePhysical;

class G4CMPInterfaceTableGenerator {
public:
  using Mat3 = G4CMPAnisotropicInterfaceSolver::Mat3;

  struct Stats {
    std::size_t totalStates = 0;
    std::size_t incidentFacingStates = 0;
    std::size_t populatedStates = 0;
    std::size_t initialSolverFailures = 0;
    std::size_t regularizedStates = 0;
    std::size_t failedStates = 0;

    G4double maximumEnergyClosure = 0.0;
    G4double maximumConditionNumber = 0.0;
  };

  G4CMPInterfaceTableGenerator(G4LatticePhysical* latticeA,
                               const Mat3& solidAToInterface,
                               G4LatticePhysical* latticeB,
                               const Mat3& solidBToInterface,
                               G4double referenceSpeed);

  G4CMPInterfaceTable Generate(const std::string& volumeAName,
                               const std::string& volumeBName,
                               const G4ThreeVector& normalAtoB, G4int nTheta,
                               G4int nPhi, Stats* stats = nullptr) const;

  static G4String GetInstalledTableDirectory();

  static G4String MakeInstalledTablePath(
      const G4String& materialA, const G4String& materialB,
      const G4RotationMatrix& solidAToInterface,
      const G4RotationMatrix& solidBToInterface);

  static G4bool TableFileExists(const G4String& filename);

private:
  G4CMPAnisotropicInterfaceSolver fSolver;

  G4LatticePhysical* fLatticeA = nullptr;
  G4LatticePhysical* fLatticeB = nullptr;
};

#endif /* G4CMPInterfaceTableGenerator_hh */