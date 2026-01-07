/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPParticleChangeForQPDiffusion.hh
/// \brief Definition of the G4CMPParticleChangeForQPDiffusion class
//
//  20260107  G4CMP-563: Delete copy operations for G4 v11 compatibility.

#ifndef G4CMPParticleChangeForQPDiffusion_h
#define G4CMPParticleChangeForQPDiffusion_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4VParticleChange.hh"

class G4DynamicParticle;

class G4CMPParticleChangeForQPDiffusion: public G4VParticleChange
{ 
public:
  // default constructor
  G4CMPParticleChangeForQPDiffusion();
  
  // destructor
  virtual ~G4CMPParticleChangeForQPDiffusion();
  
protected:
  // hide copy constructor and assignment operaor as protected
  G4CMPParticleChangeForQPDiffusion(const G4CMPParticleChangeForQPDiffusion &) = delete;
  G4CMPParticleChangeForQPDiffusion& operator=(const G4CMPParticleChangeForQPDiffusion &) = delete;


public: // with description

  // ----------------------------------------------------
  // --- the following methods are for updating G4Step -----
  // Return the pointer to the G4Step after updating the Step information
  // by using final state information of the track given by a physics
  // process
  virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
  virtual G4Step* UpdateStepForPostStep(G4Step* Step);
  // A physics process gives the final state of the particle
  // based on information of G4Track (or equivalently the PreStepPoint)
  
  virtual void Initialize(const G4Track&);
  // Initialize all propoerties by using G4Track information
  
  // ----------------------------------------------------
  //--- methods to keep information of the final state--
  //  IMPORTANT NOTE: Although the name of the class and methods are
  //   "Change", what it stores (and returns in get) are the "FINAL"
  //   values of the Position, Momentum, etc.
  
  void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
  void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);
  const G4ThreeVector* GetMomentumDirection() const;
  const G4ThreeVector* GetProposedMomentumDirection() const;
  void SetProposedMomentumDirection(const G4ThreeVector& Pfinal);
  // Get/Set theMomentumDirectionChange vector: it is the final momentum direction.
  
  const G4ThreeVector* GetPosition() const;
  void  ProposePosition(const G4ThreeVector& finalPosition);
  const G4ThreeVector* GetProposedPosition() const;
  void  SetProposedPosition(const G4ThreeVector& finalPosition);
  //  Get/Set the final position of the current particle.
  
  void ProposeVelocity(const G4double& Vfinal);
  const G4double* GetVelocity() const;
  const G4double* GetProposedVelocity() const;
  void SetProposedVelocity(const G4double& Vfinal);

  //  Get/Propose the final global/local Time
  // NOTE: DO NOT INVOKE both methods in a step
  //       Each method affects both local and global time  
  void ProposeGlobalTime(G4double t);
  void ProposeLocalTime(G4double t);

  G4double GetGlobalTime(G4double timeDelay=0.0) const;
  G4double GetLocalTime(G4double timeDelay=0.0) const;
  
public:
  virtual void DumpInfo() const;
  // for Debug
  virtual G4bool CheckIt(const G4Track&);
  
private:
  G4ThreeVector theMomentumDirection;
  //  It is the vector containing the final momentum direction
  //  after the invoked process. The application of the change
  //  of the momentum direction of the particle is not Done here.
  //  The responsibility to apply the change is up the entity
  //  which invoked the process.
  
  G4ThreeVector thePosition;
  //  The changed (final) position of a given particle.
  G4double theVelocity;
  // The effective velocity of the particle given the time step
  G4double theGlobalTime0;
  //  The global time at Initial.
  G4double theLocalTime0;
  //  The local time at Initial.
  G4double theTimeChange;
  //  The change of local time of a given particle.
};

#include "G4CMPParticleChangeForQPDiffusion.icc"
#endif

