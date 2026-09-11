/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/// \file library/include/G4CMPAnisotropicInterfaceSolver.hh
/// \brief Anisotropic elastic-interface solver.
///
/// Following the derivation in:
///
///   S. I. Rokhlin, T. K. Bolland, and L. Adler,
///   "Reflection and refraction of elastic waves on a plane interface
///   between two generally anisotropic media,"
///   J. Acoust. Soc. Am. 79, 906-918 (1986).
///
/// 20260901 C. Stone-Whitehead -- first implementation

#ifndef G4CMPAnisotropicInterfaceSolver_hh
#define G4CMPAnisotropicInterfaceSolver_hh 1

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include <array>
#include <complex>
#include <memory>
#include <vector>

class G4CMPPhononKinematics;
class G4LatticeLogical;
class G4LatticePhysical;

class G4CMPAnisotropicInterfaceSolver {
public:
  static constexpr G4int kSpaceDimension = 3;
  static constexpr G4int kNumRoots = 6;

  using Complex = std::complex<G4double>;

  struct Vec3 {
    std::array<G4double, kSpaceDimension> value{{0.0, 0.0, 0.0}};

    Vec3() = default;
    Vec3(G4double x, G4double y, G4double z) : value{{x, y, z}} {}

    G4double& operator()(G4int i) { return value[i]; }
    const G4double& operator()(G4int i) const { return value[i]; }

    G4double NormSquared() const;
    G4double Norm() const;
    G4bool Normalize();
    Vec3 Normalized() const;
  };

  struct Mat3 {
    std::array<std::array<G4double, kSpaceDimension>, kSpaceDimension> value{};

    static Mat3 Zero();
    static Mat3 Identity();

    G4double& operator()(G4int row, G4int column) { return value[row][column]; }

    const G4double& operator()(G4int row, G4int column) const {
      return value[row][column];
    }

    Mat3 Transpose() const;
  };

  enum class Side { Near, Far };

  enum SolverStatus {
    kSuccess = 0,
    kInvalidArgument,
    kInvalidIncidentState,
    kRootMatrixSingular,
    kRootEigensolverFailed,
    kRootMultiplicityFailed,
    kRootReconstructionFailed,
    kRootSelectionFailed,
    kBoundaryMatrixSingular,
    kModeIdentificationFailed,
    kGroupVelocityFailed,
    kEnergyClosureFailed
  };

  struct Outcome {
    Side side = Side::Near;
    G4int mode = -1;
    G4ThreeVector kDirectionInterface;
    G4ThreeVector groupVelocityDirectionInterface;
    G4double flux = 0.0;
    G4double probability = 0.0;
  };

  struct Result {
    G4bool valid = false;
    SolverStatus status = kInvalidArgument;
    G4String message;

    Side incidentSide = Side::Near;
    G4int incidentMode = -1;
    G4ThreeVector incidentKDirectionInterface;
    G4ThreeVector incidentSlownessInterface;
    G4double incidentFlux = 0.0;

    G4double matrixCondition = 0.0;
    G4double boundaryResidual = 0.0;
    G4double probabilitySum = 0.0;
    G4double energyClosure = 0.0;

    std::vector<Outcome> outcomes;
  };

  G4CMPAnisotropicInterfaceSolver(G4LatticePhysical* nearLattice,
                                  const Mat3& nearSolidToInterface,
                                  G4LatticePhysical* farLattice,
                                  const Mat3& farSolidToInterface,
                                  G4double scalingSpeed);

  ~G4CMPAnisotropicInterfaceSolver();

  Result Solve(Side incidentSide, G4int incidentMode,
               const G4ThreeVector& incidentKDirectionInterface) const;

  G4bool IsIncidentValid(
      Side incidentSide, G4int incidentMode,
      const G4ThreeVector& incidentKDirectionInterface) const;

  G4ThreeVector GetGroupVelocityDirection(
      Side side, G4int mode, const G4ThreeVector& kDirectionInterface) const;

  void SetRootDegeneracyTolerance(G4double value) {
    rootDegeneracyTolerance = value;
  }

  void SetRootImagTolerance(G4double value) { rootImagTolerance = value; }

  void SetFluxToleranceFraction(G4double value) {
    fluxToleranceFraction = value;
  }

  void SetClosureTolerance(G4double value) { closureTolerance = value; }

  void SetBoundaryResidualTolerance(G4double value) {
    boundaryResidualTolerance = value;
  }

private:
  struct CVec3 {
    std::array<Complex, kSpaceDimension> value{
        {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)}};

    static CVec3 Zero();

    Complex& operator()(G4int i) { return value[i]; }
    const Complex& operator()(G4int i) const { return value[i]; }

    G4double NormSquared() const;
    G4double Norm() const;
    G4bool Normalize();
  };

  struct CVec6 {
    std::array<Complex, kNumRoots> value{};

    static CVec6 Zero();

    Complex& operator()(G4int i) { return value[i]; }
    const Complex& operator()(G4int i) const { return value[i]; }
  };

  struct CMat3 {
    std::array<std::array<Complex, kSpaceDimension>, kSpaceDimension> value{};

    static CMat3 Zero();

    Complex& operator()(G4int row, G4int column) { return value[row][column]; }

    const Complex& operator()(G4int row, G4int column) const {
      return value[row][column];
    }
  };

  struct CMat6 {
    std::array<std::array<Complex, kNumRoots>, kNumRoots> value{};

    static CMat6 Zero();

    Complex& operator()(G4int row, G4int column) { return value[row][column]; }

    const Complex& operator()(G4int row, G4int column) const {
      return value[row][column];
    }
  };

  struct Root {
    Complex q = Complex(0.0, 0.0);
    CVec3 slowness;
    CVec3 polarization;
    CVec3 traction;
    G4bool propagating = false;
    G4double flux = 0.0;
  };

  struct Medium {
    G4LatticePhysical* physical = nullptr;
    const G4LatticeLogical* logical = nullptr;

    G4double density = 0.0;

    Mat3 solidToInterface = Mat3::Identity();
    Mat3 interfaceToSolid = Mat3::Identity();
    Mat3 crystalToInterface = Mat3::Identity();

    std::unique_ptr<G4CMPPhononKinematics> kinematics;
    std::array<G4double, kSpaceDimension * kSpaceDimension * kSpaceDimension *
                             kSpaceDimension>
        elasticity{};
  };

  struct LU6 {
    CMat6 lu;
    std::array<G4int, kNumRoots> rowSwap{{0, 1, 2, 3, 4, 5}};
    G4bool valid = false;
  };

  static G4int GetTensorIndex(G4int i, G4int j, G4int k, G4int l);

  static Vec3 ConvertToInternal(const G4ThreeVector& vector);
  static G4ThreeVector ConvertToG4(const Vec3& vector);

  static Vec3 Add(const Vec3& a, const Vec3& b);
  static Vec3 Scale(const Vec3& vector, G4double scale);
  static Vec3 Multiply(const Mat3& matrix, const Vec3& vector);

  static CVec3 Add(const CVec3& a, const CVec3& b);
  static CVec3 Subtract(const CVec3& a, const CVec3& b);
  static CVec3 Scale(const CVec3& vector, Complex scale);
  static Complex DotConjugate(const CVec3& a, const CVec3& b);

  Mat3 BuildCrystalToInterface(G4LatticePhysical* lattice,
                               const Mat3& solidToInterface) const;

  Medium BuildMedium(G4LatticePhysical* lattice,
                     const Mat3& solidToInterface) const;

  G4double GetElasticity(const Medium& medium, G4int i, G4int j, G4int k,
                         G4int l) const;

  CVec3 ComputeTraction(const Medium& medium, const CVec3& slowness,
                        const CVec3& polarization) const;

  G4double ComputeFlux(const Root& root) const;

  G4bool MakeBulkRoot(const Medium& medium, G4int mode,
                      const Vec3& kDirectionInterface, Root& root) const;

  void IssueWarning(const G4String& origin, const G4String& code,
                    const G4String& warningMessage) const;

  static std::array<Complex, 7> BuildSexticCoefficients(const Mat3& A0,
                                                        const Mat3& A1,
                                                        const Mat3& A2);

  G4bool SolveSextic(const std::array<Complex, 7>& coefficients,
                     std::array<Complex, kNumRoots>& roots) const;

  static CMat3 BuildQMatrix(const Mat3& A0, const Mat3& A1, const Mat3& A2,
                            Complex r);

  G4bool FindNullspace(const CMat3& matrix, G4int requestedNullity,
                       std::vector<CVec3>& basis) const;

  G4bool SolveRoots(const Medium& medium, const Vec3& tangentialSlowness,
                    std::array<Root, kNumRoots>& roots, SolverStatus& status,
                    G4String& resultMessage) const;

  G4bool SelectOutgoingRoots(const Medium& medium, Side side,
                             const std::array<Root, kNumRoots>& roots,
                             std::vector<Root>& selected,
                             G4String& resultMessage) const;

  G4bool GetKDirection(const Root& root, Vec3& kDirection) const;

  G4bool CalculateGroupVelocityDirection(const Medium& medium, G4int mode,
                                         const Vec3& kDirectionInterface,
                                         Vec3& groupVelocityDirection) const;

  G4int IdentifyMode(const Medium& medium, const Root& root) const;

  G4bool Factorize6x6(const CMat6& matrix, LU6& factorization) const;

  static CVec6 SolveFactorized6x6(const LU6& factorization,
                                  const CVec6& rightHandSide);

  G4double EstimateConditionInfinity(const CMat6& matrix,
                                     const LU6& factorization) const;

  Medium nearMedium;
  Medium farMedium;

  G4double referenceSpeed = 0.0;

  const Vec3 interfaceNormal{0.0, 0.0, 1.0};

  G4double rootDegeneracyTolerance = 1.0e-7;
  G4double rootImagTolerance = 1.0e-10;
  G4double fluxToleranceFraction = 1.0e-10;
  G4double closureTolerance = 1.0e-6;
  G4double polynomialTolerance = 5.0e-13;
  G4double nullspaceTolerance = 1.0e-9;
  G4double linearSolveTolerance = 1.0e-13;
  G4double boundaryResidualTolerance = 1.0e-11;
};

#endif /* G4CMPAnisotropicInterfaceSolver_hh */