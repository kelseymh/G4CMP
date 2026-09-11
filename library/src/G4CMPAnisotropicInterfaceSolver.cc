/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPAnisotropicInterfaceSolver.cc
/// \brief Implementation of the G4CMP anisotropic interface solver.
///
/// The continuum formulation follows:
///
///   S. I. Rokhlin, T. K. Bolland, and L. Adler,
///   "Reflection and refraction of elastic waves on a plane interface
///   between two generally anisotropic media,"
///   J. Acoust. Soc. Am. 79, 906-918 (1986).
///
/// 09012026 C. Stone-Whitehead -- first implementation
///

#include "G4CMPAnisotropicInterfaceSolver.hh"

#include "G4CMPConfigManager.hh"
#include "G4CMPPhononKinematics.hh"

#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4ios.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

// Constructor/destructor

G4CMPAnisotropicInterfaceSolver::G4CMPAnisotropicInterfaceSolver(
    G4LatticePhysical* nearLattice, const Mat3& nearSolidToInterface,
    G4LatticePhysical* farLattice, const Mat3& farSolidToInterface,
    G4double scalingSpeed)
    : referenceSpeed(scalingSpeed) {
  if (nearLattice == nullptr || farLattice == nullptr) {
    G4ExceptionDescription message;
    message << "A non-null physical lattice is required on both sides "
            << "of the interface.";
    G4Exception(
        "G4CMPAnisotropicInterfaceSolver::"
        "G4CMPAnisotropicInterfaceSolver",
        "G4CMPAniso001", FatalErrorInArgument, message);
  }

  if (!(referenceSpeed > 0.0)) {
    G4ExceptionDescription message;
    message << "The reference speed must be positive; received "
            << referenceSpeed << ".";
    G4Exception(
        "G4CMPAnisotropicInterfaceSolver::"
        "G4CMPAnisotropicInterfaceSolver",
        "G4CMPAniso002", FatalErrorInArgument, message);
  }

  // Medium helps keep rotated lattice information for the solver consistant and
  // and simplies function definitions

  nearMedium = BuildMedium(nearLattice, nearSolidToInterface);
  farMedium = BuildMedium(farLattice, farSolidToInterface);
}

G4CMPAnisotropicInterfaceSolver::~G4CMPAnisotropicInterfaceSolver() = default;

// Result helps keep track of relevant solver information for debugging
G4CMPAnisotropicInterfaceSolver::Result

G4CMPAnisotropicInterfaceSolver::Solve(
    Side incidentSide, G4int incidentMode,
    const G4ThreeVector& incidentKDirectionInterface) const {
  Result result;
  result.incidentSide = incidentSide;
  result.incidentMode = incidentMode;

  // Side allows us to encode directionality so that we can calculate
  // the physics from A->B and B->A in a convenient way
  if (incidentSide != Side::Near && incidentSide != Side::Far) {
    result.status = kInvalidArgument;
    result.message = "Invalid interface side.";

    G4ExceptionDescription message;
    message << result.message;
    G4Exception("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso007",
                FatalErrorInArgument, message);
    return result;
  }

  // Make sure incident mode corresponds to the three acoustic modes
  // ( 0, 1 or 2)
  if (incidentMode < 0 || incidentMode >= G4PhononPolarization::NUM_MODES) {
    result.status = kInvalidArgument;

    std::ostringstream stream;
    stream << "Invalid incident phonon mode " << incidentMode << ".";
    result.message = stream.str();

    G4ExceptionDescription message;
    message << result.message;
    G4Exception("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso008",
                FatalErrorInArgument, message);
    return result;
  }

  // make sure the wavevector direction is non-zero

  if (incidentKDirectionInterface.mag2() == 0.0) {
    result.status = kInvalidArgument;
    result.message = "A non-zero incident wave-vector direction is required.";

    G4ExceptionDescription message;
    message << result.message;
    G4Exception("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso009",
                FatalErrorInArgument, message);
    return result;
  }

  Vec3 incidentK = ConvertToInternal(incidentKDirectionInterface);

  if (!incidentK.Normalize()) {
    result.status = kInvalidArgument;
    result.message = "A non-zero incident wave-vector direction is required.";
    return result;
  }

  // Set the incident material information
  const Medium& incidentMedium =
      (incidentSide == Side::Near) ? nearMedium : farMedium;

  // Construct root object to keep track of incident phonon information
  Root incident;

  if (!MakeBulkRoot(incidentMedium, incidentMode, incidentK, incident)) {
    result.status = kInvalidIncidentState;
    result.message = "Failed to construct the incident bulk state.";
    return result;
  }

  result.incidentKDirectionInterface = ConvertToG4(incidentK);

  Vec3 incidentSlownessReal;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    incidentSlownessReal(i) = std::real(incident.slowness(i));
  }

  result.incidentSlownessInterface = ConvertToG4(incidentSlownessReal);

  // Confirm energy is propagating towards the interface for
  // the incident phonon

  result.incidentFlux = incident.flux;

  const G4double incidentFluxTolerance =
      fluxToleranceFraction * incidentMedium.density * referenceSpeed;

  if (incidentSide == Side::Near && !(incident.flux > incidentFluxTolerance)) {
    result.status = kInvalidIncidentState;
    result.message =
        "Near incident state does not carry energy toward the interface.";
    return result;
  }

  if (incidentSide == Side::Far && !(incident.flux < -incidentFluxTolerance)) {
    result.status = kInvalidIncidentState;
    result.message =
        "Far incident state does not carry energy toward the interface.";
    return result;
  }

  // Solve for the 6 roots that correspond to the normal component of
  // the slowness vector for the 3 reflected and 3 transmitted modes
  // for each incident direction at the interface

  Vec3 tangentialSlowness = incidentSlownessReal;
  tangentialSlowness(kSpaceDimension - 1) = 0.0;

  std::array<Root, kNumRoots> nearRoots;
  std::array<Root, kNumRoots> farRoots;

  if (!SolveRoots(nearMedium, tangentialSlowness, nearRoots, result.status,
                  result.message)) {
    return result;
  }

  if (!SolveRoots(farMedium, tangentialSlowness, farRoots, result.status,
                  result.message)) {
    return result;
  }

  std::vector<Root> nearOutgoing;
  std::vector<Root> farOutgoing;

  if (!SelectOutgoingRoots(nearMedium, Side::Near, nearRoots, nearOutgoing,
                           result.message)) {
    result.status = kRootSelectionFailed;
    return result;
  }

  if (!SelectOutgoingRoots(farMedium, Side::Far, farRoots, farOutgoing,
                           result.message)) {
    result.status = kRootSelectionFailed;
    return result;
  }

  // To get amplitudes of each root on each side, we
  // apply the displacement and traction boundary conditions

  CMat6 boundaryMatrix = CMat6::Zero();
  CVec6 rightHandSide = CVec6::Zero();

  const G4double impedanceScale =
      std::sqrt(nearMedium.density * farMedium.density) * referenceSpeed;

  for (G4int rootIndex = 0; rootIndex < kSpaceDimension; ++rootIndex) {
    for (G4int component = 0; component < kSpaceDimension; ++component) {
      boundaryMatrix(component, rootIndex) =
          nearOutgoing[rootIndex].polarization(component);

      boundaryMatrix(component + kSpaceDimension, rootIndex) =
          nearOutgoing[rootIndex].traction(component) / impedanceScale;

      boundaryMatrix(component, rootIndex + kSpaceDimension) =
          -farOutgoing[rootIndex].polarization(component);

      boundaryMatrix(component + kSpaceDimension, rootIndex + kSpaceDimension) =
          -farOutgoing[rootIndex].traction(component) / impedanceScale;
    }
  }

  const G4double rightHandSideSign = (incidentSide == Side::Near) ? -1.0 : 1.0;

  for (G4int component = 0; component < kSpaceDimension; ++component) {
    rightHandSide(component) =
        rightHandSideSign * incident.polarization(component);

    rightHandSide(component + kSpaceDimension) =
        rightHandSideSign * incident.traction(component) / impedanceScale;
  }

  // Perform LU decomposition which allows us to factorize the
  // boundary condition matrix above into an upper and lower
  // triangular matrix

  LU6 boundaryLU;

  if (!Factorize6x6(boundaryMatrix, boundaryLU)) {
    result.status = kBoundaryMatrixSingular;
    result.message =
        "RBA Eq. (30) boundary-condition matrix is numerically singular.";

    IssueWarning("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso107",
                 result.message);
    return result;
  }

  result.matrixCondition =
      EstimateConditionInfinity(boundaryMatrix, boundaryLU);

  // Solve for the displacement amplitudes
  CVec6 amplitudes = SolveFactorized6x6(boundaryLU, rightHandSide);

  // Iterative refinement

  auto multiplyBoundary = [&](const CVec6& vector) {
    CVec6 product = CVec6::Zero();

    for (G4int row = 0; row < kNumRoots; ++row) {
      for (G4int column = 0; column < kNumRoots; ++column) {
        product(row) += boundaryMatrix(row, column) * vector(column);
      }
    }

    return product;
  };

  auto vectorInfinityNorm = [&](const CVec6& vector) {
    G4double norm = 0.0;

    for (G4int i = 0; i < kNumRoots; ++i) {
      norm = std::max(norm, std::abs(vector(i)));
    }

    return norm;
  };

  G4double boundaryMatrixNorm = 0.0;

  for (G4int row = 0; row < kNumRoots; ++row) {
    G4double rowSum = 0.0;

    for (G4int column = 0; column < kNumRoots; ++column) {
      rowSum += std::abs(boundaryMatrix(row, column));
    }

    boundaryMatrixNorm = std::max(boundaryMatrixNorm, rowSum);
  }

  for (G4int refinement = 0; refinement < 4; ++refinement) {
    const CVec6 product = multiplyBoundary(amplitudes);
    CVec6 residual = CVec6::Zero();

    for (G4int i = 0; i < kNumRoots; ++i) {
      residual(i) = rightHandSide(i) - product(i);
    }

    const CVec6 correction = SolveFactorized6x6(boundaryLU, residual);

    for (G4int i = 0; i < kNumRoots; ++i) {
      amplitudes(i) += correction(i);
    }

    if (vectorInfinityNorm(correction) <=
        32.0 * std::numeric_limits<G4double>::epsilon() *
            std::max(1.0, vectorInfinityNorm(amplitudes))) {
      break;
    }
  }

  // Check that the solution satisfies the boundary conditions
  // and reject the solution if it doesn't
  {
    const CVec6 product = multiplyBoundary(amplitudes);
    CVec6 residual = CVec6::Zero();

    for (G4int i = 0; i < kNumRoots; ++i) {
      residual(i) = rightHandSide(i) - product(i);
    }

    const G4double denominator =
        boundaryMatrixNorm * vectorInfinityNorm(amplitudes) +
        vectorInfinityNorm(rightHandSide);

    result.boundaryResidual =
        vectorInfinityNorm(residual) / std::max(1.0, denominator);

    if (result.boundaryResidual > boundaryResidualTolerance) {
      result.status = kBoundaryMatrixSingular;

      std::ostringstream stream;
      stream << "RBA Eq. (30) boundary solve has excessive backward residual: "
             << result.boundaryResidual
             << ", condition_inf=" << result.matrixCondition << ".";
      result.message = stream.str();

      IssueWarning("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso107",
                   result.message);
      return result;
    }
  }
  // Calculate the probabilities from the amplitudes for each root

  const G4double incidentFluxMagnitude = std::abs(incident.flux);

  auto appendOutcomes = [&](Side side, const Medium& medium,
                            const std::vector<Root>& roots,
                            G4int amplitudeOffset) -> G4bool {
    for (G4int rootIndex = 0; rootIndex < kSpaceDimension; ++rootIndex) {
      const Root& root = roots[rootIndex];

      if (!root.propagating) {
        continue;
      }

      const G4double outwardFlux =
          (side == Side::Near) ? -root.flux : root.flux;

      if (outwardFlux <= 0.0) {
        continue;
      }
      // Convert back to quantities for G4CMP
      Outcome outcome;
      outcome.side = side;
      outcome.mode = IdentifyMode(medium, root);

      if (outcome.mode == G4PhononPolarization::UNKNOWN) {
        result.status = kModeIdentificationFailed;
        result.message =
            "Failed to identify the G4CMP mode of an outgoing root.";
        return false;
      }

      Vec3 outgoingK;

      if (!GetKDirection(root, outgoingK)) {
        result.status = kModeIdentificationFailed;
        result.message =
            "Failed to recover the phase direction of an outgoing root.";
        return false;
      }

      outcome.kDirectionInterface = ConvertToG4(outgoingK);

      Vec3 groupVelocityDirection;

      if (!CalculateGroupVelocityDirection(medium, outcome.mode, outgoingK,
                                           groupVelocityDirection)) {
        result.status = kGroupVelocityFailed;
        result.message =
            "G4CMP returned a zero outgoing group-velocity direction.";

        IssueWarning("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso108",
                     result.message);
        return false;
      }

      outcome.groupVelocityDirectionInterface =
          ConvertToG4(groupVelocityDirection);
      outcome.flux = root.flux;

      outcome.probability = std::norm(amplitudes(amplitudeOffset + rootIndex)) *
                            outwardFlux / incidentFluxMagnitude;

      result.probabilitySum += outcome.probability;
      result.outcomes.push_back(outcome);
    }

    return true;
  };

  if (!appendOutcomes(Side::Near, nearMedium, nearOutgoing, 0)) {
    return result;
  }

  if (!appendOutcomes(Side::Far, farMedium, farOutgoing, kSpaceDimension)) {
    return result;
  }

  // Ensure the probabilities sum to 1
  result.energyClosure = std::abs(1.0 - result.probabilitySum);

  if (result.energyClosure > closureTolerance) {
    result.status = kEnergyClosureFailed;

    std::ostringstream stream;
    stream << "Energy closure failed: sum=" << result.probabilitySum
           << ", closure=" << result.energyClosure
           << ", condition_inf=" << result.matrixCondition
           << ", boundary_residual=" << result.boundaryResidual << ".";
    result.message = stream.str();

    IssueWarning("G4CMPAnisotropicInterfaceSolver::Solve", "G4CMPAniso109",
                 result.message);
    return result;
  }

  result.valid = true;
  result.status = kSuccess;
  result.message = "ok";

  if (G4CMPConfigManager::GetVerboseLevel() > 0) {
    G4cout << "G4CMPAnisotropicInterfaceSolver: solved interface state with "
           << result.outcomes.size()
           << " propagating outcomes; sum(P)=" << result.probabilitySum << "."
           << G4endl;
  }

  return result;
}

// Check that the incident phonon information is valid
G4bool G4CMPAnisotropicInterfaceSolver::IsIncidentValid(
    Side incidentSide, G4int incidentMode,
    const G4ThreeVector& incidentKDirectionInterface) const {
  if (incidentSide != Side::Near && incidentSide != Side::Far) {
    return false;
  }

  if (incidentMode < 0 || incidentMode >= G4PhononPolarization::NUM_MODES) {
    return false;
  }

  if (incidentKDirectionInterface.mag2() == 0.0) {
    return false;
  }

  const Medium& medium = (incidentSide == Side::Near) ? nearMedium : farMedium;

  Root root;

  if (!MakeBulkRoot(medium, incidentMode,
                    ConvertToInternal(incidentKDirectionInterface), root)) {
    return false;
  }

  const G4double fluxTolerance =
      fluxToleranceFraction * medium.density * referenceSpeed;

  if (incidentSide == Side::Near) {
    return root.flux > fluxTolerance;
  }

  return root.flux < -fluxTolerance;
}

G4ThreeVector G4CMPAnisotropicInterfaceSolver::GetGroupVelocityDirection(
    Side side, G4int mode, const G4ThreeVector& kDirectionInterface) const {
  if (mode < 0 || mode >= G4PhononPolarization::NUM_MODES) {
    G4ExceptionDescription message;
    message << "Invalid phonon mode " << mode << ".";
    G4Exception("G4CMPAnisotropicInterfaceSolver::GetGroupVelocityDirection",
                "G4CMPAniso005", FatalErrorInArgument, message);
    return G4ThreeVector();
  }

  if (kDirectionInterface.mag2() == 0.0) {
    G4ExceptionDescription message;
    message << "A non-zero wave-vector direction is required.";
    G4Exception("G4CMPAnisotropicInterfaceSolver::GetGroupVelocityDirection",
                "G4CMPAniso006", FatalErrorInArgument, message);
    return G4ThreeVector();
  }

  const Medium& medium = (side == Side::Near) ? nearMedium : farMedium;

  Vec3 groupVelocityDirection;

  if (!CalculateGroupVelocityDirection(medium, mode,
                                       ConvertToInternal(kDirectionInterface),
                                       groupVelocityDirection)) {
    const G4String message = "G4CMP returned a zero group-velocity direction.";
    IssueWarning("G4CMPAnisotropicInterfaceSolver::GetGroupVelocityDirection",
                 "G4CMPAniso106", message);
    return G4ThreeVector();
  }

  return ConvertToG4(groupVelocityDirection);
}

// Fixed-size algebra to replace packages I previously included

G4double G4CMPAnisotropicInterfaceSolver::Vec3::NormSquared() const {
  return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

G4double G4CMPAnisotropicInterfaceSolver::Vec3::Norm() const {
  return std::sqrt(NormSquared());
}

G4bool G4CMPAnisotropicInterfaceSolver::Vec3::Normalize() {
  const G4double norm = Norm();

  if (!(norm > 0.0)) {
    return false;
  }

  for (G4double& component : value) {
    component /= norm;
  }

  return true;
}

G4CMPAnisotropicInterfaceSolver::Vec3
G4CMPAnisotropicInterfaceSolver::Vec3::Normalized() const {
  Vec3 result = *this;
  result.Normalize();
  return result;
}

G4CMPAnisotropicInterfaceSolver::Mat3
G4CMPAnisotropicInterfaceSolver::Mat3::Zero() {
  return Mat3{};
}

G4CMPAnisotropicInterfaceSolver::Mat3
G4CMPAnisotropicInterfaceSolver::Mat3::Identity() {
  Mat3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i, i) = 1.0;
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::Mat3
G4CMPAnisotropicInterfaceSolver::Mat3::Transpose() const {
  Mat3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      result(i, j) = (*this)(j, i);
    }
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::CVec3
G4CMPAnisotropicInterfaceSolver::CVec3::Zero() {
  return CVec3{};
}

G4double G4CMPAnisotropicInterfaceSolver::CVec3::NormSquared() const {
  G4double result = 0.0;

  for (const Complex& component : value) {
    result += std::norm(component);
  }

  return result;
}

G4double G4CMPAnisotropicInterfaceSolver::CVec3::Norm() const {
  return std::sqrt(NormSquared());
}

G4bool G4CMPAnisotropicInterfaceSolver::CVec3::Normalize() {
  const G4double norm = Norm();

  if (!(norm > 0.0)) {
    return false;
  }

  for (Complex& component : value) {
    component /= norm;
  }

  return true;
}

G4CMPAnisotropicInterfaceSolver::CVec6
G4CMPAnisotropicInterfaceSolver::CVec6::Zero() {
  return CVec6{};
}

G4CMPAnisotropicInterfaceSolver::CMat3
G4CMPAnisotropicInterfaceSolver::CMat3::Zero() {
  return CMat3{};
}

G4CMPAnisotropicInterfaceSolver::CMat6
G4CMPAnisotropicInterfaceSolver::CMat6::Zero() {
  return CMat6{};
}

G4int G4CMPAnisotropicInterfaceSolver::GetTensorIndex(G4int i, G4int j, G4int k,
                                                      G4int l) {
  return ((((i * kSpaceDimension) + j) * kSpaceDimension + k) *
              kSpaceDimension +
          l);
}

G4CMPAnisotropicInterfaceSolver::Vec3
G4CMPAnisotropicInterfaceSolver::ConvertToInternal(
    const G4ThreeVector& vector) {
  return Vec3(vector.x(), vector.y(), vector.z());
}

G4ThreeVector G4CMPAnisotropicInterfaceSolver::ConvertToG4(const Vec3& vector) {
  return G4ThreeVector(vector(0), vector(1), vector(2));
}

G4CMPAnisotropicInterfaceSolver::Vec3 G4CMPAnisotropicInterfaceSolver::Add(
    const Vec3& a, const Vec3& b) {
  Vec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i) = a(i) + b(i);
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::Vec3 G4CMPAnisotropicInterfaceSolver::Scale(
    const Vec3& vector, G4double scale) {
  Vec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i) = scale * vector(i);
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::Vec3 G4CMPAnisotropicInterfaceSolver::Multiply(
    const Mat3& matrix, const Vec3& vector) {
  Vec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      result(i) += matrix(i, j) * vector(j);
    }
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::CVec3 G4CMPAnisotropicInterfaceSolver::Add(
    const CVec3& a, const CVec3& b) {
  CVec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i) = a(i) + b(i);
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::CVec3
G4CMPAnisotropicInterfaceSolver::Subtract(const CVec3& a, const CVec3& b) {
  CVec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i) = a(i) - b(i);
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::CVec3 G4CMPAnisotropicInterfaceSolver::Scale(
    const CVec3& vector, Complex scale) {
  CVec3 result;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result(i) = scale * vector(i);
  }

  return result;
}

G4CMPAnisotropicInterfaceSolver::Complex
G4CMPAnisotropicInterfaceSolver::DotConjugate(const CVec3& a, const CVec3& b) {
  Complex result(0.0, 0.0);

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    result += std::conj(a(i)) * b(i);
  }

  return result;
}
G4CMPAnisotropicInterfaceSolver::Mat3
G4CMPAnisotropicInterfaceSolver::BuildCrystalToInterface(
    G4LatticePhysical* lattice, const Mat3& solidToInterface) const {
  Mat3 rotation = Mat3::Zero();

  for (G4int axisIndex = 0; axisIndex < kSpaceDimension; ++axisIndex) {
    G4ThreeVector axis(0.0, 0.0, 0.0);

    if (axisIndex == 0) {
      axis.setX(1.0);
    } else if (axisIndex == 1) {
      axis.setY(1.0);
    } else {
      axis.setZ(1.0);
    }

    lattice->RotateToSolid(axis);

    const Vec3 interfaceAxis =
        Multiply(solidToInterface, ConvertToInternal(axis));

    for (G4int row = 0; row < kSpaceDimension; ++row) {
      rotation(row, axisIndex) = interfaceAxis(row);
    }
  }

  return rotation;
}

G4CMPAnisotropicInterfaceSolver::Medium
G4CMPAnisotropicInterfaceSolver::BuildMedium(
    G4LatticePhysical* lattice, const Mat3& solidToInterface) const {
  Medium medium;

  medium.physical = lattice;
  medium.logical = lattice->GetLattice();

  if (medium.logical == nullptr) {
    G4ExceptionDescription message;
    message << "The physical lattice does not reference a logical lattice.";
    G4Exception("G4CMPAnisotropicInterfaceSolver::BuildMedium", "G4CMPAniso003",
                FatalException, message);
  }

  medium.density = medium.logical->GetDensity();

  if (!(medium.density > 0.0)) {
    G4ExceptionDescription message;
    message << "The logical-lattice density must be positive; received "
            << medium.density << ".";
    G4Exception("G4CMPAnisotropicInterfaceSolver::BuildMedium", "G4CMPAniso004",
                FatalException, message);
  }

  medium.solidToInterface = solidToInterface;
  medium.interfaceToSolid = solidToInterface.Transpose();
  medium.crystalToInterface =
      BuildCrystalToInterface(lattice, solidToInterface);

  medium.kinematics.reset(
      new G4CMPPhononKinematics(const_cast<G4LatticeLogical*>(medium.logical)));

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      for (G4int k = 0; k < kSpaceDimension; ++k) {
        for (G4int l = 0; l < kSpaceDimension; ++l) {
          G4double transformed = 0.0;

          for (G4int a = 0; a < kSpaceDimension; ++a) {
            for (G4int b = 0; b < kSpaceDimension; ++b) {
              for (G4int c = 0; c < kSpaceDimension; ++c) {
                for (G4int d = 0; d < kSpaceDimension; ++d) {
                  transformed += medium.crystalToInterface(i, a) *
                                 medium.crystalToInterface(j, b) *
                                 medium.crystalToInterface(k, c) *
                                 medium.crystalToInterface(l, d) *
                                 medium.logical->GetCijkl(a, b, c, d);
                }
              }
            }
          }

          medium.elasticity[GetTensorIndex(i, j, k, l)] = transformed;
        }
      }
    }
  }

  return medium;
}

G4double G4CMPAnisotropicInterfaceSolver::GetElasticity(const Medium& medium,
                                                        G4int i, G4int j,
                                                        G4int k,
                                                        G4int l) const {
  return medium.elasticity[GetTensorIndex(i, j, k, l)];
}

G4CMPAnisotropicInterfaceSolver::CVec3
G4CMPAnisotropicInterfaceSolver::ComputeTraction(
    const Medium& medium, const CVec3& slowness,
    const CVec3& polarization) const {
  CVec3 traction = CVec3::Zero();

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      for (G4int l = 0; l < kSpaceDimension; ++l) {
        for (G4int m = 0; m < kSpaceDimension; ++m) {
          traction(i) += GetElasticity(medium, i, j, l, m) *
                         interfaceNormal(j) * slowness(m) * polarization(l);
        }
      }
    }
  }

  return traction;
}

G4double G4CMPAnisotropicInterfaceSolver::ComputeFlux(const Root& root) const {
  return std::real(DotConjugate(root.polarization, root.traction));
}

G4bool G4CMPAnisotropicInterfaceSolver::MakeBulkRoot(
    const Medium& medium, G4int mode, const Vec3& kDirectionInterface,
    Root& root) const {
  if (mode < 0 || mode >= G4PhononPolarization::NUM_MODES) {
    return false;
  }

  if (!(kDirectionInterface.Norm() > 0.0)) {
    return false;
  }

  const Vec3 kInterface = kDirectionInterface.Normalized();
  const Vec3 kSolid = Multiply(medium.interfaceToSolid, kInterface);

  G4ThreeVector kLattice = ConvertToG4(kSolid);
  medium.physical->RotateToLattice(kLattice);

  G4ThreeVector slownessLattice =
      medium.kinematics->getSlowness(mode, kLattice);
  G4ThreeVector polarizationLattice =
      medium.kinematics->getPolarization(mode, kLattice);

  medium.physical->RotateToSolid(slownessLattice);
  medium.physical->RotateToSolid(polarizationLattice);

  const Vec3 slownessInterface =
      Multiply(medium.solidToInterface, ConvertToInternal(slownessLattice));

  Vec3 polarizationInterface =
      Multiply(medium.solidToInterface, ConvertToInternal(polarizationLattice));

  if (!polarizationInterface.Normalize()) {
    return false;
  }

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    root.slowness(i) = Complex(slownessInterface(i), 0.0);
    root.polarization(i) = Complex(polarizationInterface(i), 0.0);
  }

  root.q = root.slowness(kSpaceDimension - 1);
  root.traction = ComputeTraction(medium, root.slowness, root.polarization);
  root.propagating = true;
  root.flux = ComputeFlux(root);

  return true;
}

void G4CMPAnisotropicInterfaceSolver::IssueWarning(
    const G4String& origin, const G4String& code,
    const G4String& warningMessage) const {
  G4ExceptionDescription description;
  description << warningMessage;

  G4Exception(origin.c_str(), code.c_str(), JustWarning, description);
}

std::array<G4CMPAnisotropicInterfaceSolver::Complex, 7>
G4CMPAnisotropicInterfaceSolver::BuildSexticCoefficients(const Mat3& A0,
                                                         const Mat3& A1,
                                                         const Mat3& A2) {
  using Poly = std::array<Complex, 7>;

  auto zeroPoly = []() {
    Poly result{};
    for (Complex& coefficient : result) {
      coefficient = Complex(0.0, 0.0);
    }
    return result;
  };

  auto entry = [&](G4int row, G4int column) {
    Poly result = zeroPoly();
    result[0] = Complex(A0(row, column), 0.0);
    result[1] = Complex(A1(row, column), 0.0);
    result[2] = Complex(A2(row, column), 0.0);
    return result;
  };

  auto addPoly = [&](const Poly& a, const Poly& b) {
    Poly result = zeroPoly();
    for (G4int i = 0; i < 7; ++i) {
      result[i] = a[i] + b[i];
    }
    return result;
  };

  auto subtractPoly = [&](const Poly& a, const Poly& b) {
    Poly result = zeroPoly();
    for (G4int i = 0; i < 7; ++i) {
      result[i] = a[i] - b[i];
    }
    return result;
  };

  auto multiplyPoly = [&](const Poly& a, const Poly& b) {
    Poly result = zeroPoly();

    for (G4int i = 0; i < 7; ++i) {
      if (a[i] == Complex(0.0, 0.0)) {
        continue;
      }

      for (G4int j = 0; j + i < 7; ++j) {
        if (b[j] == Complex(0.0, 0.0)) {
          continue;
        }

        result[i + j] += a[i] * b[j];
      }
    }

    return result;
  };

  const Poly q00 = entry(0, 0);
  const Poly q01 = entry(0, 1);
  const Poly q02 = entry(0, 2);
  const Poly q10 = entry(1, 0);
  const Poly q11 = entry(1, 1);
  const Poly q12 = entry(1, 2);
  const Poly q20 = entry(2, 0);
  const Poly q21 = entry(2, 1);
  const Poly q22 = entry(2, 2);

  const Poly minor0 =
      subtractPoly(multiplyPoly(q11, q22), multiplyPoly(q12, q21));

  const Poly minor1 =
      subtractPoly(multiplyPoly(q10, q22), multiplyPoly(q12, q20));

  const Poly minor2 =
      subtractPoly(multiplyPoly(q10, q21), multiplyPoly(q11, q20));

  return addPoly(
      subtractPoly(multiplyPoly(q00, minor0), multiplyPoly(q01, minor1)),
      multiplyPoly(q02, minor2));
}

G4bool G4CMPAnisotropicInterfaceSolver::SolveSextic(
    const std::array<Complex, 7>& coefficients,
    std::array<Complex, kNumRoots>& roots) const {
  std::array<Complex, 7> normalized = coefficients;

  G4double coefficientScale = 0.0;
  for (const Complex& coefficient : normalized) {
    coefficientScale = std::max(coefficientScale, std::abs(coefficient));
  }

  if (!(coefficientScale > 0.0)) {
    return false;
  }

  for (Complex& coefficient : normalized) {
    coefficient /= coefficientScale;
  }

  if (std::abs(normalized[6]) <= polynomialTolerance) {
    return false;
  }

  const Complex leading = normalized[6];
  for (Complex& coefficient : normalized) {
    coefficient /= leading;
  }

  auto laguerre = [&](const std::vector<Complex>& polynomial, G4int degree,
                      Complex initial, G4int maximumIterations,
                      Complex& root) -> G4bool {
    static const G4double fraction[] = {0.5,  0.25, 0.75, 0.13,
                                        0.38, 0.62, 0.88, 1.0};

    Complex x = initial;

    for (G4int iteration = 0; iteration < maximumIterations; ++iteration) {
      Complex b = polynomial[degree];
      Complex d(0.0, 0.0);
      Complex f(0.0, 0.0);

      G4double errorScale = std::abs(b);
      const G4double absX = std::abs(x);

      for (G4int j = degree - 1; j >= 0; --j) {
        f = x * f + d;
        d = x * d + b;
        b = x * b + polynomial[j];
        errorScale = std::abs(b) + absX * errorScale;
      }

      if (std::abs(b) <= polynomialTolerance * std::max(1.0, errorScale)) {
        root = x;
        return true;
      }

      if (std::abs(b) == 0.0) {
        root = x;
        return true;
      }

      const Complex g = d / b;
      const Complex h = g * g - Complex(2.0, 0.0) * f / b;
      const Complex squareRoot =
          std::sqrt(Complex(static_cast<G4double>(degree - 1), 0.0) *
                    (Complex(static_cast<G4double>(degree), 0.0) * h - g * g));

      const Complex gp = g + squareRoot;
      const Complex gm = g - squareRoot;
      const Complex denominator = (std::abs(gp) > std::abs(gm)) ? gp : gm;

      Complex delta;

      if (std::abs(denominator) > 0.0) {
        delta = Complex(static_cast<G4double>(degree), 0.0) / denominator;
      } else {
        const G4double angle = 2.0 * CLHEP::pi *
                               static_cast<G4double>(iteration + 1) /
                               static_cast<G4double>(maximumIterations + 1);

        delta = std::polar(1.0 + absX, angle);
      }

      const Complex candidate = x - delta;

      if (candidate == x) {
        root = x;
        return true;
      }

      if ((iteration + 1) % 10 == 0) {
        const G4int index =
            ((iteration + 1) / 10 - 1) %
            static_cast<G4int>(sizeof(fraction) / sizeof(fraction[0]));
        x -= fraction[index] * delta;
      } else {
        x = candidate;
      }
    }

    root = x;
    return false;
  };

  std::vector<Complex> work(normalized.begin(), normalized.end());

  for (G4int degree = kNumRoots; degree >= 1; --degree) {
    Complex root;
    G4bool converged = false;

    for (G4int attempt = 0; attempt < 12 && !converged; ++attempt) {
      Complex initial(0.0, 0.0);

      if (attempt > 0) {
        G4double radius = 1.0;

        const G4double leadingMagnitude = std::abs(work[degree]);
        if (leadingMagnitude > 0.0) {
          G4double maxRatio = 0.0;

          for (G4int j = 0; j < degree; ++j) {
            maxRatio = std::max(maxRatio, std::abs(work[j]) / leadingMagnitude);
          }

          radius += maxRatio;
        }

        const G4double angle =
            2.0 * CLHEP::pi * static_cast<G4double>(attempt - 1) / 11.0 +
            0.17320508075688773;

        initial = std::polar(radius, angle);
      }

      converged = laguerre(work, degree, initial, 120, root);
    }

    if (!converged) {
      return false;
    }

    if (std::abs(std::imag(root)) <=
        rootImagTolerance * (1.0 + std::abs(root))) {
      root = Complex(std::real(root), 0.0);
    }

    roots[degree - 1] = root;

    std::vector<Complex> deflated(degree, Complex(0.0, 0.0));
    deflated[degree - 1] = work[degree];

    for (G4int j = degree - 2; j >= 0; --j) {
      deflated[j] = work[j + 1] + root * deflated[j + 1];
    }

    work = std::move(deflated);
  }

  const std::vector<Complex> original(normalized.begin(), normalized.end());

  for (Complex& root : roots) {
    Complex polished = root;

    if (laguerre(original, kNumRoots, root, 60, polished)) {
      root = polished;
    }

    if (std::abs(std::imag(root)) <=
        rootImagTolerance * (1.0 + std::abs(root))) {
      root = Complex(std::real(root), 0.0);
    }
  }

  std::sort(roots.begin(), roots.end(), [](const Complex& a, const Complex& b) {
    if (std::real(a) != std::real(b)) {
      return std::real(a) < std::real(b);
    }

    return std::imag(a) < std::imag(b);
  });

  return true;
}

G4CMPAnisotropicInterfaceSolver::CMat3

// build the complex elastic matrix where the
// determinant gives the normal-slowness roots

G4CMPAnisotropicInterfaceSolver::BuildQMatrix(const Mat3& A0, const Mat3& A1,
                                              const Mat3& A2, Complex r) {
  CMat3 result = CMat3::Zero();

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      result(i, j) = Complex(A0(i, j), 0.0) + r * Complex(A1(i, j), 0.0) +
                     r * r * Complex(A2(i, j), 0.0);
    }
  }

  return result;
}

// find allowable polarizations

G4bool G4CMPAnisotropicInterfaceSolver::FindNullspace(
    const CMat3& matrix, G4int requestedNullity,
    std::vector<CVec3>& basis) const {
  basis.clear();

  if (requestedNullity < 1 || requestedNullity > kSpaceDimension) {
    return false;
  }

  G4double matrixScale = 0.0;
  G4double imaginaryScale = 0.0;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int j = 0; j < kSpaceDimension; ++j) {
      matrixScale = std::max(matrixScale, std::abs(matrix(i, j)));
      imaginaryScale =
          std::max(imaginaryScale, std::abs(std::imag(matrix(i, j))));
    }
  }

  if (!(matrixScale > 0.0)) {
    for (G4int column = 0; column < requestedNullity; ++column) {
      CVec3 vector = CVec3::Zero();
      vector(column) = Complex(1.0, 0.0);
      basis.push_back(vector);
    }
    return true;
  }

  auto relativeResidual = [&](const CVec3& vector) -> G4double {
    CVec3 residual = CVec3::Zero();

    for (G4int i = 0; i < kSpaceDimension; ++i) {
      for (G4int j = 0; j < kSpaceDimension; ++j) {
        residual(i) += matrix(i, j) * vector(j);
      }
    }

    return residual.Norm() /
           (std::sqrt(3.0) * matrixScale * std::max(1.0, vector.Norm()));
  };

  const G4double realMatrixTolerance =
      100.0 * rootImagTolerance * std::max(1.0, matrixScale);

  if (imaginaryScale <= realMatrixTolerance) {
    G4double a[kSpaceDimension][kSpaceDimension]{};
    G4double v[kSpaceDimension][kSpaceDimension]{};

    for (G4int i = 0; i < kSpaceDimension; ++i) {
      v[i][i] = 1.0;

      for (G4int j = 0; j < kSpaceDimension; ++j) {
        a[i][j] = 0.5 * (std::real(matrix(i, j)) + std::real(matrix(j, i)));
      }
    }

    const G4double jacobiTolerance = 64.0 *
                                     std::numeric_limits<G4double>::epsilon() *
                                     std::max(1.0, matrixScale);

    for (G4int sweep = 0; sweep < 50; ++sweep) {
      G4int p = 0;
      G4int q = 1;
      G4double largest = std::abs(a[p][q]);

      for (G4int i = 0; i < kSpaceDimension; ++i) {
        for (G4int j = i + 1; j < kSpaceDimension; ++j) {
          const G4double magnitude = std::abs(a[i][j]);

          if (magnitude > largest) {
            largest = magnitude;
            p = i;
            q = j;
          }
        }
      }

      if (largest <= jacobiTolerance) {
        break;
      }

      const G4double app = a[p][p];
      const G4double aqq = a[q][q];
      const G4double apq = a[p][q];

      const G4double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
      const G4double c = std::cos(phi);
      const G4double s = std::sin(phi);

      for (G4int k = 0; k < kSpaceDimension; ++k) {
        if (k == p || k == q) {
          continue;
        }

        const G4double akp = a[k][p];
        const G4double akq = a[k][q];

        a[k][p] = a[p][k] = c * akp - s * akq;
        a[k][q] = a[q][k] = s * akp + c * akq;
      }

      a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
      a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
      a[p][q] = a[q][p] = 0.0;

      for (G4int k = 0; k < kSpaceDimension; ++k) {
        const G4double vkp = v[k][p];
        const G4double vkq = v[k][q];

        v[k][p] = c * vkp - s * vkq;
        v[k][q] = s * vkp + c * vkq;
      }
    }

    std::array<G4int, kSpaceDimension> order{{0, 1, 2}};

    std::sort(order.begin(), order.end(), [&](G4int lhs, G4int rhs) {
      return std::abs(a[lhs][lhs]) < std::abs(a[rhs][rhs]);
    });

    const G4double allowedRelativeResidual =
        std::max(1.0e-8, 100.0 * nullspaceTolerance);

    for (G4int basisIndex = 0; basisIndex < requestedNullity; ++basisIndex) {
      CVec3 vector = CVec3::Zero();
      const G4int column = order[basisIndex];

      for (G4int row = 0; row < kSpaceDimension; ++row) {
        vector(row) = Complex(v[row][column], 0.0);
      }

      if (!vector.Normalize()) {
        basis.clear();
        return false;
      }

      if (relativeResidual(vector) > allowedRelativeResidual) {
        basis.clear();
        break;
      }

      basis.push_back(vector);
    }

    if (static_cast<G4int>(basis.size()) == requestedNullity) {
      return true;
    }

    basis.clear();
  }

  if (requestedNullity == 1) {
    auto rowCross = [&](G4int first, G4int second) {
      CVec3 result = CVec3::Zero();

      result(0) = matrix(first, 1) * matrix(second, 2) -
                  matrix(first, 2) * matrix(second, 1);
      result(1) = matrix(first, 2) * matrix(second, 0) -
                  matrix(first, 0) * matrix(second, 2);
      result(2) = matrix(first, 0) * matrix(second, 1) -
                  matrix(first, 1) * matrix(second, 0);

      return result;
    };

    CVec3 best = CVec3::Zero();
    G4double bestNorm = 0.0;

    const G4int pair[3][2] = {{0, 1}, {1, 2}, {2, 0}};

    for (const auto& rows : pair) {
      CVec3 candidate = rowCross(rows[0], rows[1]);
      const G4double norm = candidate.Norm();

      if (norm > bestNorm) {
        bestNorm = norm;
        best = candidate;
      }
    }

    if (best.Normalize()) {
      const G4double allowedRelativeResidual =
          std::max(1.0e-8, 100.0 * nullspaceTolerance);

      if (relativeResidual(best) <= allowedRelativeResidual) {
        basis.push_back(best);
        return true;
      }
    }
  }

  const G4double toleranceMultiplier[] = {1.0, 10.0, 100.0};

  for (G4double multiplier : toleranceMultiplier) {
    const G4double tolerance =
        nullspaceTolerance * multiplier * std::max(1.0, matrixScale);

    std::array<std::array<Complex, kSpaceDimension>, kSpaceDimension> reduced =
        matrix.value;

    std::array<G4int, kSpaceDimension> pivotColumns{{-1, -1, -1}};
    G4int rank = 0;

    for (G4int column = 0; column < kSpaceDimension && rank < kSpaceDimension;
         ++column) {
      G4int pivotRow = -1;
      G4double pivotMagnitude = 0.0;

      for (G4int row = rank; row < kSpaceDimension; ++row) {
        const G4double magnitude = std::abs(reduced[row][column]);

        if (magnitude > pivotMagnitude) {
          pivotMagnitude = magnitude;
          pivotRow = row;
        }
      }

      if (pivotRow < 0 || pivotMagnitude <= tolerance) {
        continue;
      }

      if (pivotRow != rank) {
        std::swap(reduced[pivotRow], reduced[rank]);
      }

      const Complex pivot = reduced[rank][column];

      for (G4int j = 0; j < kSpaceDimension; ++j) {
        reduced[rank][j] /= pivot;
      }

      for (G4int row = 0; row < kSpaceDimension; ++row) {
        if (row == rank) {
          continue;
        }

        const Complex factor = reduced[row][column];

        for (G4int j = 0; j < kSpaceDimension; ++j) {
          reduced[row][j] -= factor * reduced[rank][j];
        }
      }

      pivotColumns[rank] = column;
      ++rank;
    }

    if (kSpaceDimension - rank < requestedNullity) {
      continue;
    }

    std::array<G4bool, kSpaceDimension> isPivot{{false, false, false}};

    for (G4int row = 0; row < rank; ++row) {
      isPivot[pivotColumns[row]] = true;
    }

    std::vector<CVec3> candidateBasis;

    for (G4int freeColumn = 0; freeColumn < kSpaceDimension; ++freeColumn) {
      if (isPivot[freeColumn]) {
        continue;
      }

      CVec3 vector = CVec3::Zero();
      vector(freeColumn) = Complex(1.0, 0.0);

      for (G4int row = rank - 1; row >= 0; --row) {
        const G4int pivotColumn = pivotColumns[row];
        Complex sum(0.0, 0.0);

        for (G4int column = 0; column < kSpaceDimension; ++column) {
          if (column != pivotColumn) {
            sum += reduced[row][column] * vector(column);
          }
        }

        vector(pivotColumn) = -sum;
      }

      for (const CVec3& previous : candidateBasis) {
        vector =
            Subtract(vector, Scale(previous, DotConjugate(previous, vector)));
      }

      if (vector.Normalize()) {
        candidateBasis.push_back(vector);
      }
    }

    if (static_cast<G4int>(candidateBasis.size()) < requestedNullity) {
      continue;
    }

    const G4double allowedRelativeResidual =
        std::max(1.0e-8, 100.0 * nullspaceTolerance);

    G4bool acceptable = true;

    for (G4int index = 0; index < requestedNullity; ++index) {
      if (relativeResidual(candidateBasis[index]) > allowedRelativeResidual) {
        acceptable = false;
        break;
      }
    }

    if (acceptable) {
      basis.assign(candidateBasis.begin(),
                   candidateBasis.begin() + requestedNullity);
      return true;
    }
  }

  return false;
}

G4bool G4CMPAnisotropicInterfaceSolver::SolveRoots(
    const Medium& medium, const Vec3& tangentialSlowness,
    std::array<Root, kNumRoots>& roots, SolverStatus& status,
    G4String& resultMessage) const {
  const Vec3 pbar = Scale(tangentialSlowness, referenceSpeed);

  Mat3 A0 = Mat3::Zero();
  Mat3 A1 = Mat3::Zero();
  Mat3 A2 = Mat3::Zero();

  const G4double stiffnessScale =
      medium.density * referenceSpeed * referenceSpeed;

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    for (G4int l = 0; l < kSpaceDimension; ++l) {
      for (G4int j = 0; j < kSpaceDimension; ++j) {
        for (G4int m = 0; m < kSpaceDimension; ++m) {
          const G4double cbar =
              GetElasticity(medium, i, j, l, m) / stiffnessScale;

          A2(i, l) += cbar * interfaceNormal(j) * interfaceNormal(m);

          A1(i, l) += cbar * (pbar(j) * interfaceNormal(m) +
                              interfaceNormal(j) * pbar(m));

          A0(i, l) += cbar * pbar(j) * pbar(m);
        }
      }

      if (i == l) {
        A0(i, l) -= 1.0;
      }
    }
  }

  const std::array<Complex, 7> coefficients =
      BuildSexticCoefficients(A0, A1, A2);

  G4double coefficientScale = 0.0;

  for (const Complex& coefficient : coefficients) {
    coefficientScale = std::max(coefficientScale, std::abs(coefficient));
  }

  if (!(coefficientScale > 0.0) ||
      std::abs(coefficients[6]) <= polynomialTolerance * coefficientScale) {
    status = kRootMatrixSingular;
    resultMessage =
        "The RBA sextic has a numerically vanishing leading coefficient "
        "c6=det(A2); six finite normal-slowness roots cannot be resolved.";
    IssueWarning("G4CMPAnisotropicInterfaceSolver::SolveRoots", "G4CMPAniso101",
                 resultMessage);
    return false;
  }

  std::array<Complex, kNumRoots> sexticRoots;

  if (!SolveSextic(coefficients, sexticRoots)) {
    status = kRootEigensolverFailed;
    resultMessage = "The RBA Eq. (19) sextic root solver failed to converge.";
    IssueWarning("G4CMPAnisotropicInterfaceSolver::SolveRoots", "G4CMPAniso102",
                 resultMessage);
    return false;
  }

  std::array<G4bool, kNumRoots> used{};
  G4int outputIndex = 0;

  for (G4int seed = 0; seed < kNumRoots; ++seed) {
    if (used[seed]) {
      continue;
    }

    std::vector<G4int> group;
    const Complex seedRoot = sexticRoots[seed];

    const G4double degeneracyTolerance =
        std::max(rootDegeneracyTolerance,
                 20.0 * std::sqrt(polynomialTolerance)) *
        (1.0 + std::abs(seedRoot));

    for (G4int index = seed; index < kNumRoots; ++index) {
      if (used[index]) {
        continue;
      }

      if (std::abs(sexticRoots[index] - seedRoot) <= degeneracyTolerance) {
        group.push_back(index);
        used[index] = true;
      }
    }

    const G4int multiplicity = static_cast<G4int>(group.size());

    if (multiplicity < 1 || multiplicity > kSpaceDimension) {
      status = kRootMultiplicityFailed;

      std::ostringstream stream;
      stream << "Unexpected RBA Eq. (19) root multiplicity " << multiplicity
             << ".";
      resultMessage = stream.str();

      IssueWarning("G4CMPAnisotropicInterfaceSolver::SolveRoots",
                   "G4CMPAniso103", resultMessage);
      return false;
    }

    Complex averageRoot(0.0, 0.0);

    for (G4int index : group) {
      averageRoot += sexticRoots[index];
    }

    averageRoot /= static_cast<G4double>(multiplicity);

    const G4double imaginaryTolerance =
        rootImagTolerance * (1.0 + std::abs(averageRoot));

    const G4bool propagating =
        std::abs(std::imag(averageRoot)) <= imaginaryTolerance;

    if (propagating) {
      averageRoot = Complex(std::real(averageRoot), 0.0);
    }

    const CMat3 qMatrix = BuildQMatrix(A0, A1, A2, averageRoot);

    std::vector<CVec3> nullspace;

    if (!FindNullspace(qMatrix, multiplicity, nullspace)) {
      status = kRootReconstructionFailed;
      resultMessage =
          "Failed to reconstruct the polarization nullspace for an "
          "RBA Eq. (19) normal-slowness root.";
      return false;
    }

    for (G4int basisIndex = 0; basisIndex < multiplicity; ++basisIndex) {
      if (outputIndex >= kNumRoots) {
        status = kRootReconstructionFailed;
        resultMessage =
            "Normal-slowness root reconstruction exceeded six roots.";
        IssueWarning("G4CMPAnisotropicInterfaceSolver::SolveRoots",
                     "G4CMPAniso105", resultMessage);
        return false;
      }

      Root root;
      root.q = averageRoot / referenceSpeed;

      for (G4int component = 0; component < kSpaceDimension; ++component) {
        root.slowness(component) =
            Complex(tangentialSlowness(component), 0.0) +
            root.q * Complex(interfaceNormal(component), 0.0);
      }

      root.polarization = nullspace[basisIndex];
      root.traction = ComputeTraction(medium, root.slowness, root.polarization);
      root.propagating = propagating;
      root.flux = propagating ? ComputeFlux(root) : 0.0;

      roots[outputIndex++] = root;
    }
  }

  if (outputIndex != kNumRoots) {
    status = kRootReconstructionFailed;
    resultMessage = "Failed to reconstruct all six RBA Eq. (19) roots.";
    IssueWarning("G4CMPAnisotropicInterfaceSolver::SolveRoots", "G4CMPAniso105",
                 resultMessage);
    return false;
  }

  if (G4CMPConfigManager::GetVerboseLevel() > 2) {
    G4cout << "G4CMPAnisotropicInterfaceSolver: reconstructed " << outputIndex
           << " normal-slowness roots from the RBA Eq. (19) sextic." << G4endl;
  }

  return true;
}

G4bool G4CMPAnisotropicInterfaceSolver::SelectOutgoingRoots(
    const Medium& medium, Side side, const std::array<Root, kNumRoots>& roots,
    std::vector<Root>& selected, G4String& resultMessage) const {
  selected.clear();
  selected.reserve(kSpaceDimension);

  const G4double fluxTolerance =
      fluxToleranceFraction * medium.density * referenceSpeed;

  for (const Root& root : roots) {
    if (root.propagating) {
      if (side == Side::Near && root.flux < -fluxTolerance) {
        selected.push_back(root);
      } else if (side == Side::Far && root.flux > fluxTolerance) {
        selected.push_back(root);
      }
    } else {
      const G4double imaginaryQ = std::imag(root.q);

      if (side == Side::Near && imaginaryQ < 0.0) {
        selected.push_back(root);
      } else if (side == Side::Far && imaginaryQ > 0.0) {
        selected.push_back(root);
      }
    }
  }

  if (static_cast<G4int>(selected.size()) != kSpaceDimension) {
    std::ostringstream stream;
    stream << "Expected three outgoing/decaying roots on the "
           << (side == Side::Near ? "Near" : "Far") << " side, but found "
           << selected.size()
           << ". This can occur at an exact grazing or critical state.";
    resultMessage = stream.str();
    return false;
  }

  return true;
}

G4bool G4CMPAnisotropicInterfaceSolver::GetKDirection(const Root& root,
                                                      Vec3& kDirection) const {
  if (!root.propagating) {
    return false;
  }

  for (G4int i = 0; i < kSpaceDimension; ++i) {
    kDirection(i) = std::real(root.slowness(i));
  }

  return kDirection.Normalize();
}

G4bool G4CMPAnisotropicInterfaceSolver::CalculateGroupVelocityDirection(
    const Medium& medium, G4int mode, const Vec3& kDirectionInterface,
    Vec3& groupVelocityDirection) const {
  if (!(kDirectionInterface.Norm() > 0.0)) {
    return false;
  }

  const Vec3 kInterface = kDirectionInterface.Normalized();
  const Vec3 kSolid = Multiply(medium.interfaceToSolid, kInterface);

  G4ThreeVector velocitySolid =
      medium.physical->MapKtoVDir(mode, ConvertToG4(kSolid));

  groupVelocityDirection =
      Multiply(medium.solidToInterface, ConvertToInternal(velocitySolid));

  return groupVelocityDirection.Normalize();
}

G4int G4CMPAnisotropicInterfaceSolver::IdentifyMode(const Medium& medium,
                                                    const Root& root) const {
  Vec3 kInterface;

  if (!GetKDirection(root, kInterface)) {
    return G4PhononPolarization::UNKNOWN;
  }

  const Vec3 kSolid = Multiply(medium.interfaceToSolid, kInterface);

  G4ThreeVector kLattice = ConvertToG4(kSolid);
  medium.physical->RotateToLattice(kLattice);

  G4int bestMode = G4PhononPolarization::UNKNOWN;
  G4double bestOverlap = -1.0;

  for (G4int mode = 0; mode < G4PhononPolarization::NUM_MODES; ++mode) {
    G4ThreeVector polarizationLattice =
        medium.kinematics->getPolarization(mode, kLattice);

    medium.physical->RotateToSolid(polarizationLattice);

    Vec3 polarizationInterface = Multiply(
        medium.solidToInterface, ConvertToInternal(polarizationLattice));

    if (!polarizationInterface.Normalize()) {
      continue;
    }

    CVec3 referencePolarization = CVec3::Zero();

    for (G4int i = 0; i < kSpaceDimension; ++i) {
      referencePolarization(i) = Complex(polarizationInterface(i), 0.0);
    }

    const Complex amplitude =
        DotConjugate(root.polarization, referencePolarization);

    const G4double overlap = std::norm(amplitude);

    if (overlap > bestOverlap) {
      bestOverlap = overlap;
      bestMode = mode;
    }
  }

  return bestMode;
}

G4bool G4CMPAnisotropicInterfaceSolver::Factorize6x6(const CMat6& matrix,
                                                     LU6& factorization) const {
  factorization.lu = matrix;
  factorization.rowSwap = {{0, 1, 2, 3, 4, 5}};
  factorization.valid = false;

  G4double normInfinity = 0.0;

  for (G4int row = 0; row < kNumRoots; ++row) {
    G4double rowSum = 0.0;

    for (G4int column = 0; column < kNumRoots; ++column) {
      rowSum += std::abs(matrix(row, column));
    }

    normInfinity = std::max(normInfinity, rowSum);
  }

  const G4double pivotThreshold =
      linearSolveTolerance * std::max(1.0, normInfinity);

  for (G4int column = 0; column < kNumRoots; ++column) {
    G4int pivotRow = column;
    G4double pivotMagnitude = std::abs(factorization.lu(column, column));

    for (G4int row = column + 1; row < kNumRoots; ++row) {
      const G4double magnitude = std::abs(factorization.lu(row, column));

      if (magnitude > pivotMagnitude) {
        pivotMagnitude = magnitude;
        pivotRow = row;
      }
    }

    if (pivotMagnitude <= pivotThreshold) {
      return false;
    }

    factorization.rowSwap[column] = pivotRow;

    if (pivotRow != column) {
      std::swap(factorization.lu.value[pivotRow],
                factorization.lu.value[column]);
    }

    for (G4int row = column + 1; row < kNumRoots; ++row) {
      factorization.lu(row, column) /= factorization.lu(column, column);

      const Complex multiplier = factorization.lu(row, column);

      for (G4int j = column + 1; j < kNumRoots; ++j) {
        factorization.lu(row, j) -= multiplier * factorization.lu(column, j);
      }
    }
  }

  factorization.valid = true;
  return true;
}

G4CMPAnisotropicInterfaceSolver::CVec6
G4CMPAnisotropicInterfaceSolver::SolveFactorized6x6(
    const LU6& factorization, const CVec6& rightHandSide) {
  CVec6 solution = rightHandSide;

  for (G4int column = 0; column < kNumRoots; ++column) {
    const G4int pivotRow = factorization.rowSwap[column];

    if (pivotRow != column) {
      std::swap(solution(column), solution(pivotRow));
    }
  }

  for (G4int row = 0; row < kNumRoots; ++row) {
    for (G4int column = 0; column < row; ++column) {
      solution(row) -= factorization.lu(row, column) * solution(column);
    }
  }

  for (G4int row = kNumRoots - 1; row >= 0; --row) {
    for (G4int column = row + 1; column < kNumRoots; ++column) {
      solution(row) -= factorization.lu(row, column) * solution(column);
    }

    solution(row) /= factorization.lu(row, row);
  }

  return solution;
}

G4double G4CMPAnisotropicInterfaceSolver::EstimateConditionInfinity(
    const CMat6& matrix, const LU6& factorization) const {
  G4double matrixNorm = 0.0;

  for (G4int row = 0; row < kNumRoots; ++row) {
    G4double rowSum = 0.0;

    for (G4int column = 0; column < kNumRoots; ++column) {
      rowSum += std::abs(matrix(row, column));
    }

    matrixNorm = std::max(matrixNorm, rowSum);
  }

  std::array<G4double, kNumRoots> inverseRowSums{};

  for (G4int column = 0; column < kNumRoots; ++column) {
    CVec6 unit = CVec6::Zero();
    unit(column) = Complex(1.0, 0.0);

    const CVec6 inverseColumn = SolveFactorized6x6(factorization, unit);

    for (G4int row = 0; row < kNumRoots; ++row) {
      inverseRowSums[row] += std::abs(inverseColumn(row));
    }
  }

  G4double inverseNorm = 0.0;

  for (G4double rowSum : inverseRowSums) {
    inverseNorm = std::max(inverseNorm, rowSum);
  }

  return matrixNorm * inverseNorm;
}