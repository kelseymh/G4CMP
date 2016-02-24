/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file analysis/A01/include/ChannelingDriftHit.hh
/// \brief Definition of the ChannelingDriftHit class
//
// $Id: 41fcfc8c55d9bb55cb958644e720480cb8845dd7 $
// --------------------------------------------------------------
//
#ifndef ChannelingDriftHit_h
#define ChannelingDriftHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class ChannelingDriftHit : public G4VHit
{
public:
    ChannelingDriftHit();
    ChannelingDriftHit(G4int z);
    virtual ~ChannelingDriftHit();
    ChannelingDriftHit(const ChannelingDriftHit &right);
    const ChannelingDriftHit& operator=(const ChannelingDriftHit &right);
    int operator==(const ChannelingDriftHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    inline float x();
    inline float y();
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
private:
    G4int fLayerID;
    G4double fTime;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4double fEnergy;
    
public:
    inline void SetLayerID(G4int z) { fLayerID = z; }
    inline G4int GetLayerID() const { return fLayerID; }
    inline void SetTime(G4double t) { fTime = t; }
    inline G4double GetTime() const { return fTime; }
    inline void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    inline G4ThreeVector GetLocalPos() const { return fLocalPos; }
    inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
    inline void SetEnergy(G4double energy) { fEnergy = energy; }
    inline G4double GetEnergy() const { return fEnergy; }
};

typedef G4THitsCollection<ChannelingDriftHit> ChannelingDriftHitsCollection;

extern G4Allocator<ChannelingDriftHit> ChannelingDriftHitAllocator;

inline void* ChannelingDriftHit::operator new(size_t)
{
    void* aHit;
    aHit = (void*)ChannelingDriftHitAllocator.MallocSingle();
    return aHit;
}

inline void ChannelingDriftHit::operator delete(void* aHit)
{
    ChannelingDriftHitAllocator.FreeSingle((ChannelingDriftHit*) aHit);
}

#endif


