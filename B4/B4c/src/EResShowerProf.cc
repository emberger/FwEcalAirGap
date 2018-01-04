
#include "EResShowerProf.hh"


void EResShowerProf::EResAndProfile(Double_t cut, std::string pth){

        TCanvas * ca1=new TCanvas("c2p","c2p");
        Cevent=new B4ROOTEvent();

        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
        Int_t seed = ms;
        std::unique_ptr<TRandom3> rndGen(new TRandom3(seed));
        Double_t c=cut;
        Double_t mean= 0;
        Double_t sig=0.1;



        Int_t Estart=50;
        Double_t x[11];
        Double_t eres[11];

        std::string fileName;

        for(Int_t i=0; i<10; i++) {

                fileName=pth+"/gamma"+std::to_string(Estart)+"MeV_0.4:0.4_.root";
                std::unique_ptr<TChain> ch1(new TChain("eventTree"));
                ch1->Add(fileName.c_str());
                ch1->Draw("");

                std::cout<<"hello"<<std::endl;
                std::unique_ptr<TTree> EcalTree(ch1->GetTree());
                Int_t nofEvents=EcalTree->GetEntries();



                EcalTree->SetBranchAddress("EventBranch", &Cevent);


                std::unique_ptr<TCanvas> profileC(new TCanvas("profileC", "profileC"));
                std::unique_ptr<TH1D> profile(new TH1D("showerprofile", "showerprofile",100, 0, 100));

                std::unique_ptr<TH1D> edepdistribution(new TH1D("edep","edep",100,0,Estart));


                Int_t eventcounter=0;
                Int_t coutcounter=0;
                Double_t encounter=0;



                for(Int_t j=0; j<nofEvents; j++) {
                  Double_t newGapEnergy=0;


                        EcalTree->GetEntry(j);                   //grab event from tree
                        eventcounter+=Cevent->NHits();
                        for(Int_t l=0; l<Cevent->NHits(); l++) {

                                Double_t smear=rndGen->Gaus(mean, sig);

                                Double_t smearedenergy=Cevent->Hit(l)->EnergyDeposit()+smear;
                                Cevent->Hit(l)->SetEnergyDeposit(smearedenergy);

                                if(Cevent->Hit(l)->EnergyDeposit() < c) {
                                        coutcounter++;
                                        encounter+=Cevent->Hit(l)->EnergyDeposit();

                                        Cevent->Hit(l)->SetEnergyDeposit(0.);
                                }
                                else if(Cevent->Hit(l)->EnergyDeposit()>= c) {
                                        newGapEnergy+=Cevent->Hit(l)->EnergyDeposit();
                                }
                        }

                        std::cout<<"average rejected hits: "<<coutcounter/nofEvents<<" of an average hitnumber of " <<eventcounter/nofEvents<<std::endl;
                        std::cout<<"average rejected energy "<<encounter/nofEvents<<"MeV"<<std::endl;
                        Cevent->SetGapEnergy(newGapEnergy);

                        edepdistribution->Fill(newGapEnergy);
                        Int_t cnh = Cevent->NHits();
                        Double_t integral;

                        for(Int_t i=0; i<Cevent->NumberOfInnerLayers(); i++) { //loop over all inner layers in event
                                std::string cp="InnerGapLV";

                                for(Int_t j =0; j<cnh; j++) { //loop over all hits in layer i

                                        //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                                        Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                                        if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                                profile->Fill(i, edep);

                                        } //fill histogram with hits of layer i
                                }
                                //std::cout<<photonNR<<" : "<<std::endl;


                        }
                        for(Int_t i=0; i<Cevent->NumberOfOuterLayers(); i++) { //loop over all outer layers in event
                                std::string cp="OuterGapLV";

                                for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i

                                        Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                                        if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                                profile->Fill(i+Cevent->NumberOfInnerLayers(), edep);

                                        } //fill histogram with hits of layer i
                                }

                        } //end of event

                }
                eres[i]=(edepdistribution->GetRMS(1)/edepdistribution->GetMean());

                EcalTree->GetEntry(0);
                std::unique_ptr<TLine> line_x(new TLine(Cevent->NumberOfInnerLayers()-1,0,Cevent->NumberOfInnerLayers()-1,profile->GetMaximum()));

                profileC->cd(0);
                profile->Draw();
                line_x->Draw();
                profile->GetXaxis()->SetTitle("Layerindex");
                profile->GetYaxis()->SetTitle("Energy of Run[MeV]");

                std::string path=pth+"/Showerprofile"+std::to_string(Estart)+".C";
                profileC->Print(path.c_str());
                path=pth+"/Showerprofile"+std::to_string(Estart)+".pdf";
                profileC->Print(path.c_str());



                x[i]=Estart;
                Estart+=100;
        }
        Estart=1500;

        fileName=pth+"/gamma"+std::to_string(Estart)+"MeV_0.4:0.4_.root";
        std::unique_ptr<TChain> ch1(new TChain("eventTree"));
        ch1->Add(fileName.c_str());
        ch1->Draw("");

        std::unique_ptr<TTree> EcalTree(ch1->GetTree());
        Int_t nofEvents=EcalTree->GetEntries();

        EcalTree->SetBranchAddress("EventBranch", &Cevent);

        std::unique_ptr<TCanvas> profileC(new TCanvas("profileC", "profileC"));
        std::unique_ptr<TH1D> profile(new TH1D("showerprofile", "showerprofile",100, 0, 100));

        std::unique_ptr<TH1D> edepdistribution(new TH1D("edep","edep",100,0,Estart));


        Int_t eventcounter=0;
        Int_t coutcounter=0;
        Double_t encounter=0;



        for(Int_t j=0; j<nofEvents; j++) {
                  Double_t newGapEnergy=0;
                EcalTree->GetEntry(j);             //grab event from tree
                eventcounter+=Cevent->NHits();
                for(Int_t l=0; l<Cevent->NHits(); l++) {


                        Double_t smear=rndGen->Gaus(mean, sig);

                        Double_t smearedenergy=Cevent->Hit(l)->EnergyDeposit()+smear;
                        Cevent->Hit(l)->SetEnergyDeposit(smearedenergy);


                        if(Cevent->Hit(l)->EnergyDeposit() < c) {
                                coutcounter++;
                                encounter+=Cevent->Hit(l)->EnergyDeposit();

                                Cevent->Hit(l)->SetEnergyDeposit(0.);


                        }
                        else if(Cevent->Hit(l)->EnergyDeposit()>= c) {
                                newGapEnergy+=Cevent->Hit(l)->EnergyDeposit();
                        }
                }

                std::cout<<"average rejected hits: "<<coutcounter/nofEvents<<" of an average hitnumber of " <<eventcounter/nofEvents<<std::endl;
                std::cout<<"average rejected energy "<<encounter/nofEvents<<"MeV"<<std::endl;
                Cevent->SetGapEnergy(newGapEnergy);

                edepdistribution->Fill(newGapEnergy);
                Int_t cnh = Cevent->NHits();
                Double_t integral;


                for(Int_t i=0; i<Cevent->NumberOfInnerLayers(); i++) { //loop over all inner layers in event
                        std::string cp="InnerGapLV";

                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in layer i

                                //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                                Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                                if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                        profile->Fill(i, edep);

                                } //fill histogram with hits of layer i
                        }
                        //std::cout<<photonNR<<" : "<<std::endl;


                }
                for(Int_t i=0; i<Cevent->NumberOfOuterLayers(); i++) { //loop over all outer layers in event
                        std::string cp="OuterGapLV";

                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i

                                Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                                if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                        profile->Fill(i+Cevent->NumberOfInnerLayers(), edep);

                                } //fill histogram with hits of layer i
                        }

                } //end of event

        }
        eres[10]=(edepdistribution->GetRMS(1)/edepdistribution->GetMean());

        EcalTree->GetEntry(0);
        std::unique_ptr<TLine> line_x(new TLine(Cevent->NumberOfInnerLayers()-1,0,Cevent->NumberOfInnerLayers()-1,profile->GetMaximum()));

        profileC->cd(0);
        profile->Draw();
        line_x->Draw();
        profile->GetXaxis()->SetTitle("Layerindex");
        profile->GetYaxis()->SetTitle("Energy of Run[MeV]");
        std::string path=pth+"/Showerprofile"+std::to_string(Estart)+".C";
        profileC->Print(path.c_str());
        path=pth+"/Showerprofile"+std::to_string(Estart)+".pdf";
        profileC->Print(path.c_str());

        ch1->Reset();

        x[10]=Estart;

        std::unique_ptr<TGraph> Eres(new TGraph(11, x, eres));

        TF1 * resolution = new TF1("resolution", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        ca1->cd(0);
        //Eres->Fit(resolution);
        Eres->GetXaxis()->SetTitle("Energy[MeV]");
        Eres->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
        Eres->GetYaxis()->SetTitleOffset(1.2);
        gStyle->SetOptFit(111111111);
        Eres->SetMarkerStyle(23);
        Eres->SetMarkerSize(1.1);
        Eres->Draw("ap");
        path=pth+"/EnergyResolution.C";
        ca1->Print(path.c_str());
        path=pth+"/EnergyResolution.pdf";
        ca1->Print(path.c_str());
        path=pth+"/EnergyResolution.root";
        ca1->Print(path.c_str());

}
