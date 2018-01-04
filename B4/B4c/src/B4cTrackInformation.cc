#include "B4cTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<B4cTrackInformation> aTrackInformationAllocator;

B4cTrackInformation::B4cTrackInformation()
{
        originalTrackID = 0;
        particleDefinition = 0;
        originalPosition = G4ThreeVector(0.,0.,0.);
        originalMomentum = G4ThreeVector(0.,0.,0.);
        originalEnergy = 0.;
        originalTime = 0.;
}

B4cTrackInformation::B4cTrackInformation(const G4Track* aTrack)
{
        originalTrackID = aTrack->GetTrackID();
        particleDefinition = aTrack->GetDefinition();
        originalPosition = aTrack->GetPosition();
        originalMomentum = aTrack->GetMomentum();
        originalEnergy = aTrack->GetTotalEnergy();
        originalTime = aTrack->GetGlobalTime();
}

B4cTrackInformation::B4cTrackInformation(const B4cTrackInformation* aTrackInfo)
{
        PhotonNR=aTrackInfo->PhotonNR;
        originalTrackID = aTrackInfo->originalTrackID;
        particleDefinition = aTrackInfo->particleDefinition;
        originalPosition = aTrackInfo->originalPosition;
        originalMomentum = aTrackInfo->originalMomentum;
        originalEnergy = aTrackInfo->originalEnergy;
        originalTime = aTrackInfo->originalTime;
}

B4cTrackInformation::~B4cTrackInformation(){
        ;
}

void B4cTrackInformation::Print() const
{
        G4cout
                << "Original track ID " << originalTrackID
                << " at " << originalPosition <<
                "with name: "<<particleDefinition->GetParticleName()
                <<" has momentum: "<<originalMomentum
                <<"and OriginalID: "<<PhotonNR<<G4endl;
}
