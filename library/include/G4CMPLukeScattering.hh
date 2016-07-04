/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPLukeScattering.hh
/// \brief Definition of the G4CMPLukeScattering class
//
// $Id$
//
// 20150111  New base class for both electron and hole Luke processes
// 20160110  Remerge the electron and hole subclasses into one class

#ifndef G4CMPLukeScattering_h
#define G4CMPLukeScattering_h 1

#include "globals.hh"
#include "G4CMPVDriftProcess.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include <unordered_map>

class G4VProcess;
class G4ParticleDefinition;
class G4Track;

using PDFParamTuple = std::array<G4double,3>;
using PDFDataRow = std::vector<PDFParamTuple>;
using PDFDataMatrix = std::vector<PDFDataRow>;
using PDFDataTensor = std::vector<PDFDataMatrix>;

class G4CMPLukeScattering : public G4CMPVDriftProcess {
public:
  G4CMPLukeScattering();
  virtual ~G4CMPLukeScattering();
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                        G4double prevStepSize,
                                                        G4ForceCondition* cond);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:
  // GetMeanFreePath is pure virtual, so we *have* to define it.
  virtual G4double GetMeanFreePath(const G4Track&,
                                   G4double, G4ForceCondition*) { return -1.; }

private:
  PDFDataTensor LoadDataFromLUT(const G4String& filename);
  G4double CalculateKSound(G4double mass);
  // hide assignment operator as private
  G4CMPLukeScattering(G4CMPLukeScattering&);
  G4CMPLukeScattering& operator=(const G4CMPLukeScattering& right);

  size_t ESIZE;
  size_t MACHSIZE;
  size_t THETASIZE;
  G4double EMIN;
  G4double EMAX;
  G4double MACHMIN;
  G4double MACHMAX;
  G4double THETAMIN;
  G4double THETAMAX;
  PDFDataTensor ElecGPILData;
  PDFDataTensor HoleGPILData;

#ifdef G4CMP_DEBUG
  std::ofstream output;
#endif
};

#endif	/* G4CMPLukeScattering */
