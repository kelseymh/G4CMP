/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/// \file library/src/G4CMPInterfaceTableGenerator.cc
/// \brief Anisotropic elastic-interface table generator
///
/// Uses G4CMPAnisotropicInterfaceSolver to generate a table to be used
/// during runtime
///
/// 20260901 C. Stone-Whitehead -- first implementation

#include "G4CMPInterfaceTableGenerator.hh"

#include "G4CMPConfigManager.hh"

#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace {

using InterfaceSolver = G4CMPAnisotropicInterfaceSolver;
using SolverResult = InterfaceSolver::Result;

// Because we are making a grid to use during runtime, sometimes the
// chosen directions will lie near critical directions or degeneracies.
// We can recover solutions using a small angular displacement

G4bool IsRegularizableStatus(InterfaceSolver::SolverStatus status) {
  return status == InterfaceSolver::kRootEigensolverFailed ||
         status == InterfaceSolver::kRootMultiplicityFailed ||
         status == InterfaceSolver::kRootReconstructionFailed ||
         status == InterfaceSolver::kRootSelectionFailed ||
         status == InterfaceSolver::kBoundaryMatrixSingular ||
         status == InterfaceSolver::kModeIdentificationFailed ||
         status == InterfaceSolver::kGroupVelocityFailed ||
         status == InterfaceSolver::kEnergyClosureFailed;
}

SolverResult SolveWithAngularRegularization(
    const InterfaceSolver& solver, InterfaceSolver::Side side, G4int mode,
    const G4ThreeVector& originalDirection, G4bool& recovered,
    G4double& usedPerturbation) {
  recovered = false;
  usedPerturbation = 0.0;

  SolverResult original = solver.Solve(side, mode, originalDirection);

  if (original.valid || !IsRegularizableStatus(original.status)) {
    return original;
  }

  const G4ThreeVector k = originalDirection.unit();
  const G4ThreeVector tangent1 = k.orthogonal().unit();
  const G4ThreeVector tangent2 = k.cross(tangent1).unit();

  const std::array<G4ThreeVector, 8> tangents = {
      {tangent1, -tangent1, tangent2, -tangent2, (tangent1 + tangent2).unit(),
       (tangent1 - tangent2).unit(), (-tangent1 + tangent2).unit(),
       (-tangent1 - tangent2).unit()}};

  const std::array<G4double, 6> perturbations = {
      {1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3}};

  for (const G4double epsilon : perturbations) {
    G4bool foundAtThisScale = false;
    SolverResult bestResult;
    G4double bestClosure = 1.0e300;

    for (const G4ThreeVector& tangent : tangents) {
      const G4ThreeVector candidateDirection =
          (std::cos(epsilon) * k + std::sin(epsilon) * tangent).unit();

      if (!solver.IsIncidentValid(side, mode, candidateDirection)) {
        continue;
      }

      SolverResult candidate = solver.Solve(side, mode, candidateDirection);

      if (!candidate.valid) {
        continue;
      }

      if (!foundAtThisScale || candidate.energyClosure < bestClosure) {
        foundAtThisScale = true;
        bestClosure = candidate.energyClosure;
        bestResult = candidate;
      }
    }

    if (foundAtThisScale) {
      recovered = true;
      usedPerturbation = epsilon;

      std::ostringstream message;
      message << "Regularized from solver status "
              << static_cast<G4int>(original.status)
              << " using angular perturbation " << epsilon << " rad.";

      bestResult.message = message.str();
      return bestResult;
    }
  }

  return original;
}

// prepare strings for naming reusable tables for
// computation

std::string SanitizeName(const std::string& name) {
  std::string result;

  for (std::size_t i = 0; i < name.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(name[i]);

    result += std::isalnum(c) ? static_cast<char>(c) : '_';
  }

  return result;
}

// Greatest common divisor to reduce direction vectors

long long GCD(long long a, long long b) {
  a = std::llabs(a);
  b = std::llabs(b);

  while (b != 0) {
    const long long remainder = a % b;
    a = b;
    b = remainder;
  }

  return a;
}

// use m for negative values in file names

std::string IntegerComponentLabel(long long value) {
  std::ostringstream output;

  if (value < 0) {
    output << "m" << std::llabs(value);
  } else {
    output << value;
  }

  return output.str();
}

// make direction usable for file name

std::string DirectionLabel(const G4ThreeVector& inputDirection) {
  const G4ThreeVector direction = inputDirection.unit();

  constexpr G4double scale = 1.0e6;

  long long x = static_cast<long long>(std::llround(direction.x() * scale));
  long long y = static_cast<long long>(std::llround(direction.y() * scale));
  long long z = static_cast<long long>(std::llround(direction.z() * scale));

  if (std::llabs(x) <= 1) x = 0;
  if (std::llabs(y) <= 1) y = 0;
  if (std::llabs(z) <= 1) z = 0;

  const long long divisor =
      GCD(std::llabs(x), GCD(std::llabs(y), std::llabs(z)));

  if (divisor > 0) {
    x /= divisor;
    y /= divisor;
    z /= divisor;
  }

  std::ostringstream output;
  output << IntegerComponentLabel(x) << "_" << IntegerComponentLabel(y) << "_"
         << IntegerComponentLabel(z);

  return output.str();
}

// label the file name with the orientation of the solid
// relative to the interface

std::string OrientationLabel(const G4RotationMatrix& solidToInterface) {
  // map the interface calculation back into the
  // solid frame
  const G4RotationMatrix interfaceToSolid = solidToInterface.inverse();

  // use two axes to fully describe the orientation
  const G4ThreeVector normalInSolid =
      interfaceToSolid * G4ThreeVector(0.0, 0.0, 1.0);
  const G4ThreeVector xAxisInSolid =
      interfaceToSolid * G4ThreeVector(1.0, 0.0, 0.0);

  std::ostringstream output;
  output << "N" << DirectionLabel(normalInSolid) << "_X"
         << DirectionLabel(xAxisInSolid);

  return output.str();
}

// ensure a directory exists under G4CMP install in order
// to store the table for use

void EnsureDirectoryExists(const std::string& directory) {
  struct stat info;

  if (stat(directory.c_str(), &info) == 0) {
    if (!S_ISDIR(info.st_mode)) {
      G4ExceptionDescription message;
      message << "Path exists but is not a directory: " << directory;

      G4Exception("G4CMPInterfaceTableGenerator::EnsureDirectoryExists",
                  "InterfaceTableDirectory001", FatalException, message);
    }

    return;
  }

  if (mkdir(directory.c_str(), 0755) != 0 && errno != EEXIST) {
    G4ExceptionDescription message;
    message << "Unable to create interface table directory: " << directory
            << "\nCheck that G4CMPINSTALL is writable.";

    G4Exception("G4CMPInterfaceTableGenerator::EnsureDirectoryExists",
                "InterfaceTableDirectory002", FatalException, message);
  }
}

}  // namespace

G4CMPInterfaceTableGenerator::G4CMPInterfaceTableGenerator(
    G4LatticePhysical* latticeA, const Mat3& solidAToInterface,
    G4LatticePhysical* latticeB, const Mat3& solidBToInterface,
    G4double referenceSpeed)
    : fSolver(latticeA, solidAToInterface, latticeB, solidBToInterface,
              referenceSpeed),
      fLatticeA(latticeA),
      fLatticeB(latticeB) {
  if (!fLatticeA || !fLatticeB) {
    G4ExceptionDescription message;
    message << "A non-null physical lattice is required on both sides "
            << "of the interface.";

    G4Exception(
        "G4CMPInterfaceTableGenerator::"
        "G4CMPInterfaceTableGenerator",
        "G4CMPInterfaceTable001", FatalErrorInArgument, message);
  }
}

G4CMPInterfaceTable G4CMPInterfaceTableGenerator::Generate(
    const std::string& volumeAName, const std::string& volumeBName,
    const G4ThreeVector& normalAtoB, G4int nTheta, G4int nPhi,
    Stats* statsOut) const {
  Stats stats;

  const G4int verboseLevel = G4CMPConfigManager::GetVerboseLevel();

  // Track unresolved solver failures by status

  using SolverStatus = G4CMPAnisotropicInterfaceSolver::SolverStatus;

  constexpr std::size_t kNumSolverStatuses =
      static_cast<std::size_t>(
          G4CMPAnisotropicInterfaceSolver::kEnergyClosureFailed) +
      1;

  std::array<G4long, kNumSolverStatuses> unresolvedByStatus{};
  std::array<G4long, kNumSolverStatuses> examplesPrinted{};

  auto statusName = [](SolverStatus status) -> const char* {
    switch (status) {
      case G4CMPAnisotropicInterfaceSolver::kSuccess:
        return "success";
      case G4CMPAnisotropicInterfaceSolver::kInvalidArgument:
        return "invalid argument";
      case G4CMPAnisotropicInterfaceSolver::kInvalidIncidentState:
        return "invalid incident state";
      case G4CMPAnisotropicInterfaceSolver::kRootMatrixSingular:
        return "root matrix singular";
      case G4CMPAnisotropicInterfaceSolver::kRootEigensolverFailed:
        return "sextic/root solver failed";
      case G4CMPAnisotropicInterfaceSolver::kRootMultiplicityFailed:
        return "root multiplicity failed";
      case G4CMPAnisotropicInterfaceSolver::kRootReconstructionFailed:
        return "root reconstruction failed";
      case G4CMPAnisotropicInterfaceSolver::kRootSelectionFailed:
        return "root selection failed";
      case G4CMPAnisotropicInterfaceSolver::kBoundaryMatrixSingular:
        return "boundary matrix singular";
      case G4CMPAnisotropicInterfaceSolver::kModeIdentificationFailed:
        return "mode identification failed";
      case G4CMPAnisotropicInterfaceSolver::kGroupVelocityFailed:
        return "group velocity failed";
      case G4CMPAnisotropicInterfaceSolver::kEnergyClosureFailed:
        return "energy closure failed";
    }

    return "unknown";
  };

  // Get lattice information

  const G4LatticeLogical* logicalA = fLatticeA->GetLattice();
  const G4LatticeLogical* logicalB = fLatticeB->GetLattice();

  if (!logicalA || !logicalB) {
    G4ExceptionDescription message;
    message << "A physical lattice passed to the interface table "
            << "generator does not contain a logical lattice.";

    G4Exception("G4CMPInterfaceTableGenerator::Generate",
                "G4CMPInterfaceTable002", FatalException, message);
  }

  if (normalAtoB.mag2() == 0.0) {
    G4ExceptionDescription message;
    message << "A non-zero interface normal is required.";

    G4Exception("G4CMPInterfaceTableGenerator::Generate",
                "G4CMPInterfaceTable003", FatalErrorInArgument, message);
  }

  G4CMPInterfaceTable table;

  // configure the table to be interpretable at runtime
  table.Configure(volumeAName, volumeBName, logicalA->GetName(),
                  logicalB->GetName(), normalAtoB.unit(), nTheta, nPhi);

  if (verboseLevel > 0) {
    G4cout << "\nGenerating anisotropic interface table"
           << "\n  A = " << volumeAName << " (" << logicalA->GetName() << ")"
           << "\n  B = " << volumeBName << " (" << logicalB->GetName() << ")"
           << "\n  Grid = 2 sides x 3 modes x " << nTheta << " theta x " << nPhi
           << " phi"
           << "\n  Full phase-wavevector sphere is sampled for each side/mode."
           << G4endl;
  }

  // compute for both incident directions

  for (G4int s = 0; s < 2; ++s) {
    const auto solverSide = (s == 0)
                                ? G4CMPAnisotropicInterfaceSolver::Side::Near
                                : G4CMPAnisotropicInterfaceSolver::Side::Far;

    const auto tableSide =
        (s == 0) ? G4CMPInterfaceTable::Side::A : G4CMPInterfaceTable::Side::B;

    // compute for all three possible acoustic modes

    for (G4int mode = 0; mode < G4PhononPolarization::NUM_MODES; ++mode) {
      if (verboseLevel > 1) {
        G4cout << "  side " << (s == 0 ? "A" : "B") << ", mode "
               << G4PhononPolarization::Label(mode) << " ..." << G4endl;
      }

      // sweep through incident directions

      for (G4int it = 0; it < nTheta; ++it) {
        for (G4int ip = 0; ip < nPhi; ++ip) {
          ++stats.totalStates;

          G4CMPInterfaceTable::Entry& entry = table.At(tableSide, mode, it, ip);

          entry.incidentKDir = table.GridDirection(it, ip);

          // make sure energy is propagating towards the interface

          if (!fSolver.IsIncidentValid(solverSide, mode, entry.incidentKDir)) {
            entry.status = G4CMPInterfaceTable::Status::NotIncidentFacing;
            entry.nOutcomes = 0;
            entry.probabilitySum = 0.0;
            entry.energyClosure = 0.0;
            entry.diagnostic =
                "group/energy flux does not point toward interface";
            continue;
          }

          ++stats.incidentFacingStates;

          G4bool recovered = false;
          G4double usedPerturbation = 0.0;

          // solve with the ability to recover a physical solution
          // by adjusting the angle

          auto result = SolveWithAngularRegularization(
              fSolver, solverSide, mode, entry.incidentKDir, recovered,
              usedPerturbation);

          // Recovered states are still initial failures

          if (recovered || !result.valid) {
            ++stats.initialSolverFailures;
          }

          if (recovered) {
            ++stats.regularizedStates;
          }

          // if there is no sucessful solve, store it in the table

          if (!result.valid) {
            entry.status = G4CMPInterfaceTable::Status::SolveFailed;
            entry.nOutcomes = 0;
            entry.probabilitySum = result.probabilitySum;
            entry.energyClosure = result.energyClosure;
            entry.diagnostic = result.message;

            ++stats.failedStates;

            const std::size_t statusIndex =
                static_cast<std::size_t>(result.status);

            if (statusIndex < unresolvedByStatus.size()) {
              ++unresolvedByStatus[statusIndex];

              // Print the first three failures of each class to
              // provide diagnostics output
              if (verboseLevel > 1 && examplesPrinted[statusIndex] < 3) {
                G4cout << "\nUnresolved solver failure example"
                       << "\n  status            = "
                       << statusName(result.status)
                       << "\n  status code       = "
                       << static_cast<G4int>(result.status)
                       << "\n  side              = " << (s == 0 ? "A" : "B")
                       << "\n  mode              = "
                       << G4PhononPolarization::Label(mode)
                       << "\n  theta index       = " << it << " / " << nTheta
                       << "\n  phi index         = " << ip << " / " << nPhi
                       << "\n  k direction       = " << entry.incidentKDir
                       << "\n  probability sum   = " << result.probabilitySum
                       << "\n  energy closure    = " << result.energyClosure
                       << "\n  condition_inf     = " << result.matrixCondition
                       << "\n  boundary residual = " << result.boundaryResidual
                       << "\n  message           = " << result.message
                       << G4endl;

                ++examplesPrinted[statusIndex];
              }
            } else {
              // This indicates SolverStatus changed without a corresponding
              // diagnostic stored

              G4ExceptionDescription message;
              message << "Solver returned out-of-range status code "
                      << static_cast<G4int>(result.status) << " for side "
                      << (s == 0 ? "A" : "B") << ", mode "
                      << G4PhononPolarization::Label(mode) << ", theta index "
                      << it << ", phi index " << ip << ".";

              G4Exception("G4CMPInterfaceTableGenerator::Generate",
                          "G4CMPInterfaceTable004", JustWarning, message);
            }

            continue;
          }

          // There are at most three propagating outgoing modes on each side

          if (result.outcomes.size() > 6) {
            G4ExceptionDescription message;
            message << "The anisotropic interface solver returned "
                    << result.outcomes.size()
                    << " propagating outcomes; at most six are allowed.";

            G4Exception("G4CMPInterfaceTableGenerator::Generate",
                        "G4CMPInterfaceTable005", FatalException, message);
          }

          // get result into a format for use with the table

          entry.status = G4CMPInterfaceTable::Status::Populated;
          entry.nOutcomes = static_cast<G4int>(result.outcomes.size());
          entry.probabilitySum = result.probabilitySum;
          entry.energyClosure = result.energyClosure;
          entry.diagnostic = result.message;

          if (recovered) {
            std::ostringstream diagnostic;
            diagnostic
                << result.message
                << " Original grid direction retained; solution evaluated "
                << usedPerturbation << " rad away.";

            entry.diagnostic = diagnostic.str();
          }

          // Put the solutions into the table
          for (G4int io = 0; io < entry.nOutcomes; ++io) {
            const auto& source = result.outcomes[static_cast<std::size_t>(io)];
            auto& destination = entry.outcomes[static_cast<std::size_t>(io)];

            destination.outgoingSide =
                (source.side == G4CMPAnisotropicInterfaceSolver::Side::Near)
                    ? G4CMPInterfaceTable::Side::A
                    : G4CMPInterfaceTable::Side::B;
            destination.mode = source.mode;
            destination.kDir = source.kDirectionInterface;
            destination.vgDir = source.groupVelocityDirectionInterface;
            destination.probability = source.probability;
          }

          stats.maximumEnergyClosure =
              std::max(stats.maximumEnergyClosure, result.energyClosure);
          stats.maximumConditionNumber =
              std::max(stats.maximumConditionNumber, result.matrixCondition);

          ++stats.populatedStates;
        }
      }
    }
  }

  if (verboseLevel > 0) {
    G4cout << "\nGeneration statistics:"
           << "\n  total grid states       = " << stats.totalStates
           << "\n  incident-facing states  = " << stats.incidentFacingStates
           << "\n  populated states        = " << stats.populatedStates
           << "\n  initial solver failures = " << stats.initialSolverFailures
           << "\n  regularized states      = " << stats.regularizedStates
           << "\n  unresolved failures     = " << stats.failedStates
           << "\n  max energy closure      = " << stats.maximumEnergyClosure
           << "\n  max condition number    = " << stats.maximumConditionNumber
           << G4endl;
  }

  if (verboseLevel > 0) {
    G4cout << "\nUnresolved failure breakdown:" << G4endl;
  }

  G4long histogramTotal = 0;

  for (std::size_t i = 0; i < unresolvedByStatus.size(); ++i) {
    if (unresolvedByStatus[i] == 0) {
      continue;
    }

    const auto status = static_cast<SolverStatus>(i);

    if (verboseLevel > 0) {
      G4cout << "  " << statusName(status) << " = " << unresolvedByStatus[i]
             << G4endl;
    }

    histogramTotal += unresolvedByStatus[i];
  }

  if (verboseLevel > 0) {
    G4cout << "  histogram total = " << histogramTotal
           << "\n  stats total     = " << stats.failedStates << G4endl;
  }

  if (histogramTotal != stats.failedStates) {
    G4ExceptionDescription message;
    message << "Unresolved-failure histogram total (" << histogramTotal
            << ") does not match stats.failedStates (" << stats.failedStates
            << ").";

    G4Exception("G4CMPInterfaceTableGenerator::Generate",
                "G4CMPInterfaceTable006", JustWarning, message);
  }

  if (statsOut) {
    *statsOut = stats;
  }

  return table;
}

// helpers to install the table into the InterfaceTables folder next to
// CrystalMaps

G4String G4CMPInterfaceTableGenerator::GetInstalledTableDirectory() {
  const char* install = std::getenv("G4CMPINSTALL");

  if (!install || !*install) {
    G4Exception("G4CMPInterfaceTableGenerator::GetInstalledTableDirectory",
                "InterfaceTableDirectory003", FatalException,
                "G4CMPINSTALL is not defined.");

    return "";
  }

  const std::string directory = std::string(install) + "/InterfaceTables";

  EnsureDirectoryExists(directory);

  return directory;
}

G4String G4CMPInterfaceTableGenerator::MakeInstalledTablePath(
    const G4String& materialA, const G4String& materialB,
    const G4RotationMatrix& solidAToInterface,
    const G4RotationMatrix& solidBToInterface) {
  const std::string directory = GetInstalledTableDirectory();

  const std::string nameA = SanitizeName(materialA);
  const std::string nameB = SanitizeName(materialB);

  const std::string orientationA = OrientationLabel(solidAToInterface);
  const std::string orientationB = OrientationLabel(solidBToInterface);

  std::ostringstream filename;
  filename << directory << "/" << nameA << "_" << orientationA << "__" << nameB
           << "_" << orientationB << "_InterfaceTable.dat";

  return filename.str();
}

G4bool G4CMPInterfaceTableGenerator::TableFileExists(const G4String& filename) {
  struct stat info;

  if (stat(filename.c_str(), &info) != 0) {
    return false;
  }

  return S_ISREG(info.st_mode);
}