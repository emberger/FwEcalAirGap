
#ifndef B4cTrackInformation_h
#define B4cTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class B4cTrackInformation : public G4VUserTrackInformation
{
public:
B4cTrackInformation();
B4cTrackInformation(const G4Track* aTrack);
B4cTrackInformation(const B4cTrackInformation* aTrackInfo);
virtual ~B4cTrackInformation();

inline void *operator new(size_t);
inline void operator delete(void *aTrackInfo);
inline int operator ==(const B4cTrackInformation& right) const
{
        return (this==&right);
}

void Print() const;

private:
G4int PhotonNR;
G4int originalTrackID;
G4ParticleDefinition* particleDefinition;
G4ThreeVector originalPosition;
G4ThreeVector originalMomentum;
G4double originalEnergy;
G4double originalTime;

public:


inline G4int SetOriginalPhotonNumber(G4int pn){     //set and get number of direct pion decay secondaries
        PhotonNR=pn;
}
inline G4int GetOriginalPhotonNumber(){
        return PhotonNR;
}




inline G4int GetOriginalTrackID() const {
        return originalTrackID;
}
inline G4ParticleDefinition* GetOriginalParticle() const {
        return particleDefinition;
}
inline G4ThreeVector GetOriginalPosition() const {
        return originalPosition;
}
inline G4ThreeVector GetOriginalMomentum() const {
        return originalMomentum;
}
inline G4double GetOriginalEnergy() const {
        return originalEnergy;
}
inline G4double GetOriginalTime() const {
        return originalTime;
}
};

extern G4Allocator<B4cTrackInformation> aTrackInformationAllocator;

inline void* B4cTrackInformation::operator new(size_t)
{
        void* aTrackInfo;
        aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
        return aTrackInfo;
}

inline void B4cTrackInformation::operator delete(void *aTrackInfo)
{
        aTrackInformationAllocator.FreeSingle((B4cTrackInformation*)aTrackInfo);
}

#endif
