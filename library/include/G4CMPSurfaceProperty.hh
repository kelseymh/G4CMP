/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20160831  M. Kelsey -- Add optional electrode geometry class
// 20170525  M. Kelsey -- Add default "rule of five" copy/move operators
// 20170627  M. Kelsey -- Return non-const pointers, for functional use
// 20190806  M. Kelsey -- Add local data for frequency-dependent scattering
//		probabilities, and computation functions.
// 20200601  G4CMP-206: Need thread-local copies of electrode pointers

#ifndef G4CMPSurfaceProperty_h
#define G4CMPSurfaceProperty_h 1

#include "G4SurfaceProperty.hh"
#include "G4MaterialPropertiesTable.hh"
#include <vector>
#include <map>


class G4CMPVElectrodePattern;


class G4CMPSurfaceProperty : public G4SurfaceProperty {
public:
  // Empty constructor. Users must call at least one of the FillPropertiesTable
  // member functions. But, really, you shouldn't use this. It's dangerous and
  // I don't know why I put it here at all.
  G4CMPSurfaceProperty(const G4String& name,
                       G4SurfaceType stype = dielectric_dielectric);

  //Full constructor
  G4CMPSurfaceProperty(const G4String& name,
                       G4double qAbsProb, // Prob. to absorb charge carrier
                       G4double qReflProb, // If not absorbed, prob to reflect
                       G4double eMinK, //Min wave number to absorb electron
                       G4double hMinK, //Min wave number to absorb hole
                       G4double pAbsProb, // Prob. to absorb phonon
                       G4double pReflProb, // If not absorbed, prob to reflect
                       G4double pSpecProb, //Prob. of specular reflection
                       G4double pMinK, //Min wave number to absorb phonon
                       G4SurfaceType stype = dielectric_dielectric);

  virtual ~G4CMPSurfaceProperty();

  G4CMPSurfaceProperty(const G4CMPSurfaceProperty&) = default;
  G4CMPSurfaceProperty(G4CMPSurfaceProperty&&) = default;
  G4CMPSurfaceProperty& operator=(const G4CMPSurfaceProperty&) = default;
  G4CMPSurfaceProperty& operator=(G4CMPSurfaceProperty&&) = default;

  G4bool operator==(const G4SurfaceProperty &right) const;
  G4bool operator!=(const G4SurfaceProperty &right) const;

  // Accessors for charge-pair and phonon boundary parameters
  // NOTE:  Must return non-const pointer as Tables can't be const
  G4MaterialPropertiesTable* GetChargeMaterialPropertiesTablePointer() const {
    return const_cast<G4MaterialPropertiesTable*>(&theChargeMatPropTable);
  }

  G4MaterialPropertiesTable* GetPhononMaterialPropertiesTablePointer() const {
    return const_cast<G4MaterialPropertiesTable*>(&thePhononMatPropTable);
  }

  // NOTE:  These return by value because Tables can't be const
  G4MaterialPropertiesTable
  GetChargeMaterialPropertiesTable() const { return theChargeMatPropTable; }

  G4MaterialPropertiesTable
  GetPhononMaterialPropertiesTable() const { return thePhononMatPropTable; }

  // Accessors to fill charge-pair and phonon boundary parameters
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable *mpt);
  void SetChargeMaterialPropertiesTable(G4MaterialPropertiesTable& mpt);
  void SetPhononMaterialPropertiesTable(G4MaterialPropertiesTable& mpt);

  void FillChargeMaterialPropertiesTable(G4double qAbsProb, G4double qReflProb,
                                         G4double eMinK,    G4double hMinK);

  void FillPhononMaterialPropertiesTable(G4double pAbsProb,  G4double pReflProb,
                                         G4double pSpecProb, G4double pMinK);

  // Accessors to fill phonon surface interaction parametrizations
  void AddSurfaceAnharmonicCutoff(G4double freqMax) { anharmonicMaxFreq = freqMax; }
  void AddSurfaceDiffuseCutoff(G4double freqDiff) { diffuseMaxFreq = freqDiff; }

  // For polynomial coeffients, units can be factored out and passed separately
  void AddSurfaceAnharmonicCoeffs(const std::vector<G4double>& coeff,
				  G4double freqUnits=0.) {
    SaveCoeffs(anharmonicCoeffs, coeff, freqUnits);
  }

  void AddDiffuseReflectionCoeffs(const std::vector<G4double>& coeff,
				  G4double freqUnits=0.) {
    SaveCoeffs(diffuseCoeffs, coeff, freqUnits);
  }

  void AddSpecularReflectionCoeffs(const std::vector<G4double>& coeff,
				   G4double freqUnits=0.) {
    SaveCoeffs(specularCoeffs, coeff, freqUnits);
  }

  // Functions to compute reflection probabilities vs. frequency
  G4double AnharmonicReflProb(G4double freq) const;
  G4double DiffuseReflProb(G4double freq) const;
  G4double SpecularReflProb(G4double freq) const;

  // Complex electrode geometries
  void SetChargeElectrode(G4CMPVElectrodePattern* cel);
  void SetPhononElectrode(G4CMPVElectrodePattern* pel);

  // Accessors, used by worker threads
  G4CMPVElectrodePattern* GetChargeElectrode() const;
  G4CMPVElectrodePattern* GetPhononElectrode() const;

  virtual void DumpInfo() const;	// To be implemented

  void AddScatteringProperties(G4double AnhCutoff, G4double DiffCutoff, 
	const std::vector<G4double>& AnhCoeffs, const std::vector<G4double>& DiffCoeffs,
	const std::vector<G4double>& SpecCoeffs, G4double AnhFreqUnits, G4double DiffFreqUnits,
  G4double SpecFreqUnits); 
  //Sets anharmonic, diffuse, and specular reflection properties

protected:
  // These args should be const, but G4MaterialPropertiesTables is silly.
  G4bool IsValidChargePropTable(G4MaterialPropertiesTable& propTab) const;
  G4bool IsValidPhononPropTable(G4MaterialPropertiesTable& propTab) const;

  void SaveCoeffs(std::vector<G4double>& buffer,
		  const std::vector<G4double>& coeff, G4double units);

  G4double ExpandCoeffsPoly(G4double freq, const std::vector<G4double>& coeff) const;

protected:
  G4MaterialPropertiesTable theChargeMatPropTable;
  G4MaterialPropertiesTable thePhononMatPropTable;

  G4CMPVElectrodePattern* theChargeElectrode;
  G4CMPVElectrodePattern* thePhononElectrode;

  // Frequency dependent phonon surface scattering parameters
  G4double anharmonicMaxFreq;		// Max frequency for anharmonic decay
  G4double diffuseMaxFreq;		// Limit frqeuency for diffuse scatter

  // Polynomial coefficients in frequency for each kind of scattering
  std::vector<G4double> anharmonicCoeffs;
  std::vector<G4double> diffuseCoeffs;
  std::vector<G4double> specularCoeffs;

  // These lists will be pre-allocated, with values entered by thread
  mutable std::map<G4int, G4CMPVElectrodePattern*> workerChargeElectrode;
  mutable std::map<G4int, G4CMPVElectrodePattern*> workerPhononElectrode;

  // These args should be const, but G4MaterialPropertiesTables is silly.
//   G4bool IsValidChargePropTable(G4MaterialPropertiesTable& propTab) const;
//   G4bool IsValidPhononPropTable(G4MaterialPropertiesTable& propTab) const;
};

#endif
