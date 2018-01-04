#include "B4cTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "B4cTrackInformation.hh"
#include "B4cPionManager.hh"


B4cTrackingAction::B4cTrackingAction() : G4UserTrackingAction()
{
        ;
}

B4cTrackingAction::~B4cTrackingAction()
{
        ;
}

void B4cTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{



        if(aTrack->GetTrackID()==1 && aTrack->GetDefinition()->GetParticleName()=="pi0") {
                isPion=true;
                isdecayed=false;
                std::cout<<"Pion"<<std::endl;
        }


        if(aTrack->GetTrackID()==1 && isPion==false) {
                isdecayed=true;
                //G4cout<<"TrackID2: "<<aTrack->GetTrackID()<<G4endl;
                //G4cout<<"true1"<<G4endl;
                //G4cout<<aTrack->GetVolume()->GetName()<<G4endl;

                B4cTrackInformation* anInfo = new B4cTrackInformation(aTrack);
                anInfo->SetOriginalPhotonNumber(1);
                //G4cout<<anInfo->GetOriginalPhotonNumber()<<G4endl;

                //G4cout<<"anInfo:"<<G4endl;
                //anInfo->Print();
                G4Track* theTrack = (G4Track*)aTrack;
                theTrack->SetUserInformation(anInfo);

        }





        if(aTrack->GetTrackID()==2&& isPion==true)
        {
                isdecayed=true;
                //G4cout<<"TrackID1: "<<aTrack->GetTrackID()<<G4endl;
                //G4cout<<"true1"<<G4endl;
                //G4cout<<aTrack->GetVolume()->GetName()<<G4endl;
                GetEInst().SetEnergyPh1(aTrack->GetTotalEnergy());
                GetEInst().SetMomPh1(aTrack->GetMomentum());
                B4cTrackInformation* anInfo = new B4cTrackInformation(aTrack);
                anInfo->SetOriginalPhotonNumber(1);
                G4cout<<aTrack->GetTotalEnergy()<<G4endl;
                //G4cout<<"anInfo:"<<G4endl;
                //anInfo->Print();
                G4Track* theTrack = (G4Track*)aTrack;
                theTrack->SetUserInformation(anInfo);
                //G4cout<<"true2"<<G4endl;


        }

        if(aTrack->GetTrackID()==3 && isPion==true)
        {
                isdecayed=true;
                //G4cout<<"TrackID2: "<<aTrack->GetTrackID()<<G4endl;
                //G4cout<<"true1"<<G4endl;
                //  G4cout<<aTrack->GetVolume()->GetName()<<G4endl;
                GetEInst().SetEnergyPh2(aTrack->GetTotalEnergy());
                GetEInst().SetMomPh2(aTrack->GetMomentum());
                G4cout<<aTrack->GetTotalEnergy()<<G4endl;
                B4cTrackInformation* anInfo = new B4cTrackInformation(aTrack);
                anInfo->SetOriginalPhotonNumber(2);
                //G4cout<<anInfo->GetOriginalPhotonNumber()<<G4endl;

                //G4cout<<"anInfo:"<<G4endl;
                //anInfo->Print();
                G4Track* theTrack = (G4Track*)aTrack;
                theTrack->SetUserInformation(anInfo);
                //G4cout<<"true2"<<G4endl;


        }


}

void B4cTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
        //G4cout<<"true3"<<G4endl;
        if(isdecayed) {
                G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
                if(secondaries)
                {
                        B4cTrackInformation * info;
                        //G4cout<<"true4"<<G4endl;
                        if(aTrack->GetUserInformation()) {
                                info = (B4cTrackInformation*)(aTrack->GetUserInformation());
                        }
                        else{
                                info= new B4cTrackInformation();
                                GetEInst().SetTD();

                        }
                        //G4cout<<"true4.1"<<G4endl;
                        size_t nSeco = secondaries->size();
                        //G4cout<<nSeco<<G4endl;
                        if(nSeco>0)
                        {
                                //G4cout<<"true5"<<G4endl;
                                for(size_t i=0; i<nSeco; i++)
                                {
                                        B4cTrackInformation* infoNew = new B4cTrackInformation(info);
                                        //G4cout<<"true5.1"<<G4endl;
                                        //infoNew->Print();
                                        (*secondaries)[i]->SetUserInformation(infoNew);
                                        //G4cout<<"true5.2"<<G4endl;
                                }
                        }
                }
        }
}
