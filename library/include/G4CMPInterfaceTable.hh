
/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/// \file library/include/G4CMPInterfaceTable.hh
/// \brief Anisotropic elastic-interface table
///
/// Save the table and lookup phonon information from runtime to get
/// relevant entries in the table
///
/// 20260901 C. Stone-Whitehead -- first implementation

#ifndef G4CMPInterfaceTable_hh
#define G4CMPInterfaceTable_hh 1

#include "G4ThreeVector.hh"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "globals.hh"

class G4CMPInterfaceTable {
public:
  enum class Side : G4int { A = 0, B = 1 };

  enum class Status : G4int {
    NotIncidentFacing = 0,
    Populated = 1,
    SolveFailed = 2
  };

  struct Outcome {
    Side outgoingSide = Side::A;
    G4int mode = -1;

    G4ThreeVector kDir;
    G4ThreeVector vgDir;

    G4double probability = 0.0;
  };

  struct Entry {
    Status status = Status::NotIncidentFacing;

    G4ThreeVector incidentKDir;

    std::array<Outcome, 6> outcomes;
    G4int nOutcomes = 0;

    G4double probabilitySum = 0.0;
    G4double energyClosure = 1.0;

    std::string diagnostic;
  };

  G4CMPInterfaceTable() = default;

  void Configure(const std::string& volumeA, const std::string& volumeB,
                 const std::string& latticeA, const std::string& latticeB,
                 const G4ThreeVector& normalAtoB, G4int nTheta, G4int nPhi);

  G4int GetNTheta() const { return fNTheta; }
  G4int GetNPhi() const { return fNPhi; }

  const std::string& GetVolumeA() const { return fVolumeA; }
  const std::string& GetVolumeB() const { return fVolumeB; }
  const std::string& GetLatticeA() const { return fLatticeA; }
  const std::string& GetLatticeB() const { return fLatticeB; }

  const G4ThreeVector& GetNormalAtoB() const { return fNormalAtoB; }

  std::size_t Size() const { return fEntries.size(); }

  Entry& At(Side side, G4int mode, G4int iTheta, G4int iPhi);

  const Entry& At(Side side, G4int mode, G4int iTheta, G4int iPhi) const;

  G4ThreeVector GridDirection(G4int iTheta, G4int iPhi) const;

  const Entry* LookupNearest(Side side, G4int mode, const G4ThreeVector& kDir,
                             G4int* matchedTheta = nullptr,
                             G4int* matchedPhi = nullptr) const;

  void Save(const std::string& filename) const;

  static G4CMPInterfaceTable Load(const std::string& filename);

  void WriteSummaryCSV(const std::string& filename) const;

private:
  std::size_t Index(Side side, G4int mode, G4int iTheta, G4int iPhi) const;

  void ValidateIndices(Side side, G4int mode, G4int iTheta, G4int iPhi) const;

  std::string fVolumeA;
  std::string fVolumeB;
  std::string fLatticeA;
  std::string fLatticeB;

  G4ThreeVector fNormalAtoB{0.0, 0.0, 1.0};

  G4int fNTheta = 0;
  G4int fNPhi = 0;

  std::vector<Entry> fEntries;
};

#endif /* G4CMPInterfaceTable_hh */