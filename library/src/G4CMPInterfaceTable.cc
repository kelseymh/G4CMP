/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/// \file library/src/G4CMPInterfaceTable.cc
/// \brief Anisotropic elastic-interface table
///
/// Save the table and lookup phonon information from runtime to get
/// relevant entries in the table
///
/// 20260901 C. Stone-Whitehead -- first implementation

#include "G4CMPInterfaceTable.hh"

#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4PhononPolarization.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

namespace {

// make text labels for diagnostics

const char* SideLabel(G4CMPInterfaceTable::Side side) {
  return side == G4CMPInterfaceTable::Side::A ? "A" : "B";
}

const char* StatusLabel(G4CMPInterfaceTable::Status status) {
  switch (status) {
    case G4CMPInterfaceTable::Status::NotIncidentFacing:
      return "not_incident_facing";
    case G4CMPInterfaceTable::Status::Populated:
      return "populated";
    case G4CMPInterfaceTable::Status::SolveFailed:
      return "solve_failed";
  }

  return "unknown";
}

std::string ModeLabel(G4int mode) {
  if (mode == G4PhononPolarization::Long) return "L";

  if (mode == G4PhononPolarization::TransSlow) return "ST";

  if (mode == G4PhononPolarization::TransFast) return "FT";

  return "UNKNOWN";
}
}  // namespace

// store metadata (will need to change how the volume is stored for wide use
// without regenerating), normalize normal vector, set the dimensions of the
// grid and allocate storage for the grid entries

void G4CMPInterfaceTable::Configure(const std::string& volumeA,
                                    const std::string& volumeB,
                                    const std::string& latticeA,
                                    const std::string& latticeB,
                                    const G4ThreeVector& normalAtoB,
                                    G4int nTheta, G4int nPhi) {
  if (nTheta < 2) {
    G4ExceptionDescription message;
    message << "nTheta must be greater than or equal to 2; received " << nTheta
            << ".";

    G4Exception("G4CMPInterfaceTable::Configure", "G4CMPInterfaceTable001",
                FatalErrorInArgument, message);
    return;
  }

  if (nPhi < 1) {
    G4ExceptionDescription message;
    message << "nPhi must be greater than or equal to 1; received " << nPhi
            << ".";

    G4Exception("G4CMPInterfaceTable::Configure", "G4CMPInterfaceTable002",
                FatalErrorInArgument, message);
    return;
  }

  if (normalAtoB.mag2() == 0.) {
    G4ExceptionDescription message;
    message << "A non-zero interface normal is required.";

    G4Exception("G4CMPInterfaceTable::Configure", "G4CMPInterfaceTable003",
                FatalErrorInArgument, message);
    return;
  }

  fVolumeA = volumeA;
  fVolumeB = volumeB;

  fLatticeA = latticeA;
  fLatticeB = latticeB;

  fNormalAtoB = normalAtoB.unit();

  fNTheta = nTheta;
  fNPhi = nPhi;

  const std::size_t count =
      static_cast<std::size_t>(2) *
      static_cast<std::size_t>(G4PhononPolarization::NUM_MODES) *
      static_cast<std::size_t>(fNTheta) * static_cast<std::size_t>(fNPhi);

  fEntries.assign(count, Entry{});
}

// convert necessary inputs into a flat vector index

std::size_t G4CMPInterfaceTable::Index(Side side, G4int mode, G4int iTheta,
                                       G4int iPhi) const {
  const std::size_t s = static_cast<std::size_t>(side == Side::A ? 0 : 1);

  return ((s * G4PhononPolarization::NUM_MODES +
           static_cast<std::size_t>(mode)) *
              static_cast<std::size_t>(fNTheta) +
          static_cast<std::size_t>(iTheta)) *
             static_cast<std::size_t>(fNPhi) +
         static_cast<std::size_t>(iPhi);
}

// check the logical inputs before passing
// to index

void G4CMPInterfaceTable::ValidateIndices(Side side, G4int mode, G4int iTheta,
                                          G4int iPhi) const {
  if (side != Side::A && side != Side::B) {
    G4ExceptionDescription message;
    message << "Invalid interface side.";

    G4Exception("G4CMPInterfaceTable::ValidateIndices",
                "G4CMPInterfaceTable004", FatalErrorInArgument, message);
    return;
  }

  if (mode < 0 || mode >= G4PhononPolarization::NUM_MODES) {
    G4ExceptionDescription message;
    message << "Invalid phonon mode " << mode << ".";

    G4Exception("G4CMPInterfaceTable::ValidateIndices",
                "G4CMPInterfaceTable005", FatalErrorInArgument, message);
    return;
  }

  if (iTheta < 0 || iTheta >= fNTheta) {
    G4ExceptionDescription message;
    message << "Theta index " << iTheta << " is outside [0," << fNTheta - 1
            << "].";

    G4Exception("G4CMPInterfaceTable::ValidateIndices",
                "G4CMPInterfaceTable006", FatalErrorInArgument, message);
    return;
  }

  if (iPhi < 0 || iPhi >= fNPhi) {
    G4ExceptionDescription message;
    message << "Phi index " << iPhi << " is outside [0," << fNPhi - 1 << "].";

    G4Exception("G4CMPInterfaceTable::ValidateIndices",
                "G4CMPInterfaceTable007", FatalErrorInArgument, message);
    return;
  }
}

// return a mutable table entry after validation

G4CMPInterfaceTable::Entry& G4CMPInterfaceTable::At(Side side, G4int mode,
                                                    G4int iTheta, G4int iPhi) {
  ValidateIndices(side, mode, iTheta, iPhi);

  return fEntries[Index(side, mode, iTheta, iPhi)];
}

// just for inspection

const G4CMPInterfaceTable::Entry& G4CMPInterfaceTable::At(Side side, G4int mode,
                                                          G4int iTheta,
                                                          G4int iPhi) const {
  ValidateIndices(side, mode, iTheta, iPhi);

  return fEntries[Index(side, mode, iTheta, iPhi)];
}

// return the wavevector direction at a specified point

G4ThreeVector G4CMPInterfaceTable::GridDirection(G4int iTheta,
                                                 G4int iPhi) const {
  if (iTheta < 0 || iTheta >= fNTheta || iPhi < 0 || iPhi >= fNPhi) {
    G4ExceptionDescription message;
    message << "Invalid angular grid indices: theta=" << iTheta
            << ", phi=" << iPhi << ".";

    G4Exception("G4CMPInterfaceTable::GridDirection", "G4CMPInterfaceTable008",
                FatalErrorInArgument, message);
    return G4ThreeVector();
  }

  const G4double theta = CLHEP::pi * static_cast<G4double>(iTheta) /
                         static_cast<G4double>(fNTheta - 1);

  const G4double phi =
      CLHEP::twopi * static_cast<G4double>(iPhi) / static_cast<G4double>(fNPhi);

  return G4ThreeVector(std::sin(theta) * std::cos(phi),
                       std::sin(theta) * std::sin(phi), std::cos(theta));
}

// lookup the nearest entry in the table
// to that of the incoming phonon

const G4CMPInterfaceTable::Entry* G4CMPInterfaceTable::LookupNearest(
    Side side, G4int mode, const G4ThreeVector& inputK, G4int* matchedTheta,
    G4int* matchedPhi) const {
  if (fNTheta < 2 || fNPhi < 1) {
    return nullptr;
  }

  if (inputK.mag2() == 0.) {
    return nullptr;
  }

  const G4ThreeVector k = inputK.unit();

  const G4double z = std::max(-1.0, std::min(1.0, static_cast<double>(k.z())));

  const G4double theta = std::acos(z);

  G4double phi = std::atan2(k.y(), k.x());

  if (phi < 0.) phi += CLHEP::twopi;

  G4int iTheta = static_cast<G4int>(
      std::llround(theta / CLHEP::pi * static_cast<G4double>(fNTheta - 1)));

  G4int iPhi = static_cast<G4int>(
      std::llround(phi / CLHEP::twopi * static_cast<G4double>(fNPhi)));

  iTheta = std::max(0, std::min(fNTheta - 1, iTheta));

  iPhi %= fNPhi;

  if (iPhi < 0) iPhi += fNPhi;

  const Entry& nearest = At(side, mode, iTheta, iPhi);

  if (nearest.status == Status::Populated) {
    if (matchedTheta) *matchedTheta = iTheta;

    if (matchedPhi) *matchedPhi = iPhi;

    return &nearest;
  }

  // if we are very near grazing incidence or a failed entry, search a
  // small area around the angle to ensure it is still incident facing
  // or to find a close solution

  const Entry* best = nullptr;
  G4double bestDot = -2.;
  G4int bestTheta = -1;
  G4int bestPhi = -1;

  constexpr G4int maxRadius = 4;

  for (G4int radius = 1; radius <= maxRadius; ++radius) {
    for (G4int dt = -radius; dt <= radius; ++dt) {
      const G4int ti = iTheta + dt;

      if (ti < 0 || ti >= fNTheta) continue;

      for (G4int dp = -radius; dp <= radius; ++dp) {
        if (std::abs(dt) != radius && std::abs(dp) != radius) continue;

        G4int pi = (iPhi + dp) % fNPhi;

        if (pi < 0) pi += fNPhi;

        const Entry& candidate = At(side, mode, ti, pi);

        if (candidate.status != Status::Populated) continue;

        const G4double dot = k.dot(candidate.incidentKDir.unit());

        if (dot > bestDot) {
          bestDot = dot;
          best = &candidate;
          bestTheta = ti;
          bestPhi = pi;
        }
      }
    }

    if (best) break;
  }

  if (best) {
    if (matchedTheta) *matchedTheta = bestTheta;

    if (matchedPhi) *matchedPhi = bestPhi;
  }

  return best;
}

// save the table to the specified filepath

void G4CMPInterfaceTable::Save(const std::string& filename) const {
  std::ofstream out(filename);

  if (!out) {
    G4ExceptionDescription message;
    message << "Unable to open interface table file for writing: " << filename;

    G4Exception("G4CMPInterfaceTable::Save", "G4CMPInterfaceTable009",
                FatalException, message);
    return;
  }

  out << std::setprecision(17);

  out << "G4CMP_ANISO_INTERFACE_TABLE 1\n";

  out << "volumeA " << std::quoted(fVolumeA) << "\n";

  out << "volumeB " << std::quoted(fVolumeB) << "\n";

  out << "latticeA " << std::quoted(fLatticeA) << "\n";

  out << "latticeB " << std::quoted(fLatticeB) << "\n";

  out << "normal " << fNormalAtoB.x() << " " << fNormalAtoB.y() << " "
      << fNormalAtoB.z() << "\n";

  out << "ntheta " << fNTheta << "\n";

  out << "nphi " << fNPhi << "\n";

  out << "BEGIN_ENTRIES\n";

  for (G4int s = 0; s < 2; ++s) {
    const Side side = s == 0 ? Side::A : Side::B;

    for (G4int mode = 0; mode < G4PhononPolarization::NUM_MODES; ++mode) {
      for (G4int it = 0; it < fNTheta; ++it) {
        for (G4int ip = 0; ip < fNPhi; ++ip) {
          const Entry& e = At(side, mode, it, ip);

          out << "E " << s << " " << mode << " " << it << " " << ip << " "
              << static_cast<G4int>(e.status) << " " << e.nOutcomes << " "
              << e.incidentKDir.x() << " " << e.incidentKDir.y() << " "
              << e.incidentKDir.z() << " " << e.probabilitySum << " "
              << e.energyClosure << " " << std::quoted(e.diagnostic) << "\n";

          for (G4int io = 0; io < e.nOutcomes; ++io) {
            const Outcome& o = e.outcomes[static_cast<std::size_t>(io)];

            out << "O " << static_cast<G4int>(o.outgoingSide) << " " << o.mode
                << " " << o.probability << " " << o.kDir.x() << " "
                << o.kDir.y() << " " << o.kDir.z() << " " << o.vgDir.x() << " "
                << o.vgDir.y() << " " << o.vgDir.z() << "\n";
          }
        }
      }
    }
  }

  out << "END_ENTRIES\n";
}

// load the table for use at runtime

G4CMPInterfaceTable G4CMPInterfaceTable::Load(const std::string& filename) {
  std::ifstream in(filename);

  if (!in) {
    G4ExceptionDescription message;
    message << "Unable to open interface table file for reading: " << filename;

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable010",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  std::string magic;
  G4int version = 0;

  in >> magic >> version;

  if (magic != "G4CMP_ANISO_INTERFACE_TABLE" || version != 1) {
    G4ExceptionDescription message;
    message << "Unsupported interface table format in file " << filename << ".";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable011",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  std::string key;

  std::string volumeA;
  std::string volumeB;
  std::string latticeA;
  std::string latticeB;

  G4ThreeVector normal;

  G4int nTheta = 0;
  G4int nPhi = 0;

  in >> key >> std::quoted(volumeA);

  if (key != "volumeA") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected volumeA.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable012",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  in >> key >> std::quoted(volumeB);

  if (key != "volumeB") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected volumeB.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable013",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  in >> key >> std::quoted(latticeA);

  if (key != "latticeA") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected latticeA.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable014",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  in >> key >> std::quoted(latticeB);

  if (key != "latticeB") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected latticeB.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable015",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;

  in >> key >> nx >> ny >> nz;

  if (key != "normal") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected normal.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable016",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  normal = G4ThreeVector(nx, ny, nz);

  in >> key >> nTheta;

  if (key != "ntheta") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected ntheta.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable017",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  in >> key >> nPhi;

  if (key != "nphi") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected nphi.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable018",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  in >> key;

  if (key != "BEGIN_ENTRIES") {
    G4ExceptionDescription message;
    message << "Malformed interface table: expected BEGIN_ENTRIES.";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable019",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  G4CMPInterfaceTable table;

  table.Configure(volumeA, volumeB, latticeA, latticeB, normal, nTheta, nPhi);

  const std::size_t expected = table.Size();

  std::size_t entriesRead = 0;

  while (in >> key) {
    if (key == "END_ENTRIES") break;

    if (key != "E") {
      G4ExceptionDescription message;
      message << "Malformed interface table: expected E or END_ENTRIES.";

      G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable020",
                  FatalException, message);
      return G4CMPInterfaceTable();
    }

    G4int s = 0;
    G4int mode = 0;
    G4int it = 0;
    G4int ip = 0;
    G4int status = 0;
    G4int nOut = 0;

    G4double kx = 0.;
    G4double ky = 0.;
    G4double kz = 0.;

    G4double pSum = 0.;
    G4double closure = 1.;

    std::string diagnostic;

    in >> s >> mode >> it >> ip >> status >> nOut >> kx >> ky >> kz >> pSum >>
        closure >> std::quoted(diagnostic);

    if (!in) {
      G4ExceptionDescription message;
      message << "Malformed interface table E record.";

      G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable021",
                  FatalException, message);
      return G4CMPInterfaceTable();
    }

    if (nOut < 0 || nOut > 6) {
      G4ExceptionDescription message;
      message << "Invalid outcome count " << nOut
              << " in interface table E record.";

      G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable022",
                  FatalException, message);
      return G4CMPInterfaceTable();
    }

    const Side side = s == 0 ? Side::A : Side::B;

    Entry& e = table.At(side, mode, it, ip);

    e.status = static_cast<Status>(status);

    e.nOutcomes = nOut;

    e.incidentKDir = G4ThreeVector(kx, ky, kz);

    e.probabilitySum = pSum;

    e.energyClosure = closure;

    e.diagnostic = diagnostic;

    for (G4int io = 0; io < nOut; ++io) {
      in >> key;

      if (key != "O") {
        G4ExceptionDescription message;
        message << "Malformed interface table: expected O record.";

        G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable023",
                    FatalException, message);
        return G4CMPInterfaceTable();
      }

      G4int outSide = 0;
      G4int outMode = -1;

      G4double p = 0.;

      G4double okx = 0.;
      G4double oky = 0.;
      G4double okz = 0.;

      G4double vgx = 0.;
      G4double vgy = 0.;
      G4double vgz = 0.;

      in >> outSide >> outMode >> p >> okx >> oky >> okz >> vgx >> vgy >> vgz;

      if (!in) {
        G4ExceptionDescription message;
        message << "Malformed interface table O record.";

        G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable024",
                    FatalException, message);
        return G4CMPInterfaceTable();
      }

      Outcome& o = e.outcomes[static_cast<std::size_t>(io)];

      o.outgoingSide = outSide == 0 ? Side::A : Side::B;

      o.mode = outMode;

      o.probability = p;

      o.kDir = G4ThreeVector(okx, oky, okz);

      o.vgDir = G4ThreeVector(vgx, vgy, vgz);
    }

    ++entriesRead;
  }

  if (entriesRead != expected) {
    G4ExceptionDescription message;
    message << "Interface table contains " << entriesRead
            << " entries; expected " << expected << ".";

    G4Exception("G4CMPInterfaceTable::Load", "G4CMPInterfaceTable025",
                FatalException, message);
    return G4CMPInterfaceTable();
  }

  return table;
}

// write diagnostic summary

void G4CMPInterfaceTable::WriteSummaryCSV(const std::string& filename) const {
  std::ofstream out(filename);

  if (!out) {
    G4ExceptionDescription message;
    message << "Unable to open summary CSV file for writing: " << filename;

    G4Exception("G4CMPInterfaceTable::WriteSummaryCSV",
                "G4CMPInterfaceTable026", FatalException, message);
    return;
  }

  out << std::setprecision(17);

  out << "side,mode,itheta,iphi,"
      << "status,nout,"
      << "kx,ky,kz,"
      << "probability_sum,energy_closure\n";

  for (G4int s = 0; s < 2; ++s) {
    const Side side = s == 0 ? Side::A : Side::B;

    for (G4int mode = 0; mode < G4PhononPolarization::NUM_MODES; ++mode) {
      for (G4int it = 0; it < fNTheta; ++it) {
        for (G4int ip = 0; ip < fNPhi; ++ip) {
          const Entry& e = At(side, mode, it, ip);

          out << SideLabel(side) << "," << ModeLabel(mode) << "," << it << ","
              << ip << "," << StatusLabel(e.status) << "," << e.nOutcomes << ","
              << e.incidentKDir.x() << "," << e.incidentKDir.y() << ","
              << e.incidentKDir.z() << "," << e.probabilitySum << ","
              << e.energyClosure << "\n";
        }
      }
    }
  }
}
