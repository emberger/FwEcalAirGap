#include "TROOTAnalysis.hh"

TROOTAnalysis::TROOTAnalysis(){
}


TROOTAnalysis::TROOTAnalysis(std::unique_ptr<TChain> &ch,Double_t prodist) :

        EcalTree(ch->GetTree()),

        D(new TCanvas("DC","Direction vector difference")),
        dx1(new TH1D("DeltaX_photon1","DeltaX_photon1", 1000,-3,3)),
        dy1(new TH1D("DeltaY_photon1","DeltaY_photon1", 1000,-3,3)),
        dz1(new TH1D("DeltaZ_photon1","DeltaZ_photon1", 1000,-3,3)),
        dx2(new TH1D("DeltaX_photon1_noPCA","DeltaX_photon1_noPCA", 1000,-3,3)),
        dy2(new TH1D("DeltaY_photon1_noPCA","DeltaY_photon1_noPCA", 1000,-3,3)),
        dz2(new TH1D("DeltaZ_photon1_noPCA","DeltaZ_photon1_noPCA", 1000,-3,3)),
        projectionCODR(new TCanvas("ProjectionCODR", "ProjectionCODR")),
        projectionCPCA(new TCanvas("ProjectionCPCA", "ProjectionCPCA")),
        correlationC(new TCanvas("CorrelationC", "CorrelationC")),

        angledifferenceCPCA(new TCanvas("angledifferenceCPCA","angledifferenceCPCA")),
        angledifferenceCODR(new TCanvas("angledifferenceCODR", "angledifferenceCODR")),

        angledifference_pca(new TH1D("AngleDiffPCA", "AngleDiffPCA", 1000,0,5)),
        angledifference_minimization(new TH1D("AngleDiffODR", "AngleDiffODR",1000,0,5))

{

        std::cout << "assigned Tree" << '\n';

        nofEntries=EcalTree->GetEntries();

        Cevent=new B4ROOTEvent();
        EcalTree->SetBranchAddress("EventBranch", &Cevent);

        EcalTree->GetEntry(0);

        InnerGapThickness=Cevent->InnerGapThickness();
        InnerAbsoThickness=Cevent->InnerAbsoThickness();
        InnertiledimX=Cevent->InnerTilesizeX();
        InnertiledimY=Cevent->InnerTilesizeY();
        nofInnerLayers=Cevent->NumberOfInnerLayers();
        Innerabsofirst=Cevent->GetInnerAbsFirst();


        OuterGapThickness=Cevent->OuterGapThickness();
        OuterAbsoThickness=Cevent->OuterAbsoThickness();
        OutertiledimX=Cevent->OuterTilesizeX();
        OutertiledimY=Cevent->OuterTilesizeY();
        nofOuterLayers=Cevent->NumberOfOuterLayers();
        Outerabsofirst=Cevent->GetOuterAbsFirst();

        gunposition=Cevent->GunPos();

        calsizeXY=Cevent->calsizeXY();
        PVesselThickness=Cevent->PVesselThickness();


        true_direction = Cevent->MomentumPh1();
        true_energy = Cevent->EnergyPrimary();

        NofRejectedHits=0;

        // std::cout<<"Absorber thickness is "<<AbsoThickness<<"mm"<<std::endl;
        // std::cout<<"Scintillator thickness is "<<GapThickness<<"mm"<<std::endl;
        // std::cout<<"Tilesize is "<<tiledimX<<"*"<<tiledimY<<"mm^2"<<std::endl;
        // std::cout << std::boolalpha;
        // std::cout<<"AbsorberFirst: "<<absofirst<<std::endl;

        if(InnertiledimX<OutertiledimX) {
                histsizeX=calsizeXY/InnertiledimX;
                histsizeY=calsizeXY/InnertiledimY;
                histsizeZ=nofInnerLayers+nofOuterLayers;
        }
        else{
                histsizeX=calsizeXY/OutertiledimX;
                histsizeY=calsizeXY/OutertiledimY;
                histsizeZ=nofInnerLayers+nofOuterLayers;
        }
        std::unique_ptr<TH2D> _h1(new TH2D("h1", "h1", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5));
        std::unique_ptr<TH2D> _h2(new TH2D("h2", "h2", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5));

        h1=std::move(_h1);
        h2=std::move(_h2);

        Double_t k=(-prodist-gunposition.Z())/true_direction.Z();

        TVector3 Intersect(gunposition.X()+k*true_direction.X(),
                           gunposition.Y()+k*true_direction.Y(),
                           gunposition.Z()+k*true_direction.Z());
        projection_true=Intersect;

        std::unique_ptr<TH2D> _projection_pca(new TH2D("Projection_pca", "Projection_pca",
                                                       2000,Intersect.X()-1000.0,Intersect.X()+1000.0,2000,Intersect.Y()-1000,Intersect.Y()+1000));
        projection_pca=std::move(_projection_pca);

        std::unique_ptr<TH2D> _projection_minimization(new TH2D("Projection_minimization", "Projection_minimization",
                                                                2000,Intersect.X()-1000.0,Intersect.X()+1000.0,2000,Intersect.Y()-1000.0,Intersect.Y()+1000.0));
        projection_minimization=std::move(_projection_minimization);

        std::unique_ptr<TH2D> _projection_corrleationX(new TH2D("projection_correlationX", "projection_correlationX",
                                                                2000,Intersect.X()-500.0,Intersect.X()+500.0,2000,Intersect.Y()-500.0,Intersect.Y()+500.0));
        projection_correlationX=std::move(_projection_corrleationX);

        std::unique_ptr<TH2D> _projection_corrleationY(new TH2D("projection_correlationY", "projection_correlationY",
                                                                2000,Intersect.X()-500.0,Intersect.X()+500.0,2000,Intersect.Y()-500.0,Intersect.Y()+500.0));
        projection_correlationY=std::move(_projection_corrleationY);

        std::unique_ptr<TH2D> _projection_correlationDeltaR(new TH2D("projection_correlationDeltaR","projection_correlationDeltaR", 401,-100,300,401,-100,300));
        projection_correlationDeltaR=std::move(_projection_correlationDeltaR);

        std::unique_ptr<TH1D> _projection_correlationDeltaRdiff(new TH1D("projection_DeltaRdifference","projection_DealtaRdiff", 1001,-500,500));
        projection_correlationDeltaRdiff=std::move(_projection_correlationDeltaRdiff);


        std::unique_ptr<TLine> _line_x(new TLine(Intersect.X(),Intersect.Y()-1000,Intersect.X(), Intersect.Y()));
        std::unique_ptr<TLine> _line_y(new TLine(Intersect.X()-1000.0,Intersect.Y(),Intersect.X(), Intersect.Y()));

        line_x=std::move(_line_x);
        line_y=std::move(_line_y);

        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
        Int_t seed = ms;
        std::unique_ptr<TRandom3> r(new TRandom3(seed));
        rndGen=std::move(r);

        std::cout<<"Start Analysis"<<std::endl;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

TROOTAnalysis::~TROOTAnalysis(){
}


Int_t TROOTAnalysis::GetNofEntries(){
        return nofEntries;
}

void TROOTAnalysis::ApplySmearingAndCut(Double_t c, Double_t mean, Double_t sig){

        for(Int_t i; i<Cevent->NHits(); i++) {
                Double_t smear=rndGen->Gaus(mean, sig);
                Double_t smearedenergy=Cevent->Hit(i)->EnergyDeposit()+smear;
                Cevent->Hit(i)->SetEnergyDeposit(smearedenergy);


                if(Cevent->Hit(i)->EnergyDeposit() < c) {
                        Cevent->Hit(i)->SetEnergyDeposit(0.);
                }

        }

}


void TROOTAnalysis::PrintEdep(){

        // for(Int_t i=0; i<Cevent->NHits(); i++) {
        //         std::cout<<Cevent->Hit(i)->EnergyDeposit()<<std::endl;
        // }
}

void TROOTAnalysis::CalcAngularResolution(std::string pth, std::string title){

        TCanvas * ca1=new TCanvas("c1","c1");
        TCanvas * ca2=new TCanvas("c2","c2");


        std::string fileName;
        Int_t Estart=50;
        Double_t X[11];
        Double_t Y[11];

        for(Int_t i = 0; i<10; i++) {

                fileName=pth+"/Projection_E:"+std::to_string(Estart)+"MeV_.root";
                std::cout<<fileName<<std::endl;
                TFile * f1=new TFile(fileName.c_str());

                TCanvas * c1 =(TCanvas*)f1->Get("ProjectionC");
                TH2D * hi1 = (TH2D*)c1->GetPrimitive("Projection_minimization");

                TCanvas dummyC;

                Double_t rmsx=hi1->GetRMS(1);
                //Double_t rmsy=h1->GetRMS(2);
                X[i]=Estart;
                Y[i]=rmsx/1150;
                Estart+=100;

        }

        Estart=1500;
        fileName=pth+"/Projection_E:"+std::to_string(Estart)+"MeV_.root";

        TFile * f1=new TFile(fileName.c_str());
        TCanvas * c1 =(TCanvas*)f1->Get("ProjectionC");
        TH2D * hi1 = (TH2D*)c1->GetPrimitive("Projection_minimization");

        TCanvas dummyC;
        Double_t rmsx=hi1->GetRMS(1);
        //Double_t rmsy=h1->GetRMS(2);

        X[10]=Estart;
        Y[10]=rmsx/1150;

        TGraph * g1 = new TGraph(11, X, Y);

        Int_t len=0;


        g1->SetTitle(title.c_str());
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");
        //TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");




        //g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1->GetXaxis()->SetTitle("Energy[MeV]");
        g1->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");

        gStyle->SetOptFit();

        ca1->cd(0);
        g1->Draw("A*");

        std::string path=pth+"/"+title+".root";
        ca1->Print(path.c_str());

        path=pth+"/"+title+".C";
        ca1->Print(path.c_str());

        path=pth+"/"+title+".pdf";
        ca1->Print(path.c_str());


}

//void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event
// TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");
//
// if(nofEntries<pev+1) {
//         std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
//         return;
// }
//
//
//
// //for(Int_t i=0;i<nent; i++){
//
// EcalTree->GetEntry(pev);
// Int_t pnh=Cevent->NHits();
//
// for(Int_t j=0; j<pnh; j++) {
//         h->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
// }
//
// h->GetXaxis()->SetTitle("X");
// h->GetYaxis()->SetTitle("Y");
// h->GetZaxis()->SetTitle("Z");
//
// plotcanvas1->cd();
// h->Draw("BOX");
//
// std::string Path="/home/iwsatlas1/emberger/FuerFrank/PhotonEventDisplay/Event";
// std::string nr = std::to_string(pev);
//
// std::string extension = ".C";
// std::string PlotPath=Path+nr+extension;
// plotcanvas1->Print(PlotPath.c_str());
//
// extension = ".pdf";
// PlotPath=Path+nr+extension;
// plotcanvas1->Print(PlotPath.c_str());
// h->Reset();
//}

// void TROOTAnalysis::plotEventPion(Int_t pev){     //plot 3DHisto of selected event
//         TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");



//hA= new TH3D("ECalEvent1","ECalEvent1",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
//hB= new TH3D("ECalEvent2","ECalEvent2",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
//         if(nofEntries<pev+1) {
//                 std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
//                 return;
//         }
//
//
//
//         //for(Int_t i=0;i<nent; i++){
//
//         EcalTree->GetEntry(pev);
//         Int_t pnh=Cevent->NHits();
//
//         for(Int_t j=0; j<pnh; j++) {
//                 if(Cevent->Hit(j)->PhotNr()==1) {
//                         hA->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
//                 }
//         }
//
//         //TColor * col1= new TColor(1,0,0,"red");
//
//         hA->GetXaxis()->SetTitle("X");
//         hA->GetYaxis()->SetTitle("Y");
//         hA->GetZaxis()->SetTitle("Z");
//         hA->SetMarkerColor(kRed);
//         hA->SetMarkerStyle(20);
//         hA->SetMarkerSize(0.3);
//         hA->Draw("");
//
//         for(Int_t j=0; j<pnh; j++) {
//                 if(Cevent->Hit(j)->PhotNr()==2) {
//                         hB->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
//                 }
//         }
//
//         //TColor * col2 = new TColor(0,0,1,"blue");
//         hB->SetMarkerColor(kBlue);
//         hB->SetMarkerStyle(20);
//         hB->SetMarkerSize(0.3);
//         hB->GetXaxis()->SetTitle("X");
//         hB->GetYaxis()->SetTitle("Y");
//         hB->GetZaxis()->SetTitle("Z");
//         hB->Draw("SAME");
//
//         std::string Path="/home/iwsatlas1/emberger/PionPlots/Event";
//         std::string nr = std::to_string(pev);
//
//         std::string extension = ".C";
//         std::string PlotPath=Path+nr+extension;
//         plotcanvas1->Print(PlotPath.c_str());
//
//         // extension = ".png";
//         // PlotPath=Path+nr+extension;
//         // plotcanvas1->Print(PlotPath.c_str());
//
//         extension = ".pdf";
//         PlotPath=Path+nr+extension;
//         plotcanvas1->Print(PlotPath.c_str());
//
//         hB->Reset();
//         hA->Reset();
// }


//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotRMSx(){
        TCanvas * rms=new TCanvas("rmsx_over_E");


        // Double_t y[10]={0.5298,0.3507,0.3624,0.2524,0.222,0.1861,0.228,0.1751,0.1425,0.1481}; // 0.0 slope
        // Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        //
        Double_t y[10]={0.5398,0.3428,0.3087,0.2319,0.2117,0.1933,0.2243,0.1716,0.1411,0.1588};                                                                      //0.4slope
        Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        TGraph *gRMS1=new TGraph(10, x, y );
        gRMS1->Fit(rms_x);
        gStyle->SetOptFit(11111111);
        gRMS1->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS1->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS1->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(0);
        gRMS1->Draw("A*");
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

// void TROOTAnalysis::SampleFromHIst(){
//
//         TCanvas * dummy = new TCanvas();
//         dummy->Divide(2,2,0.01,0.01);
//         TFile *f = new TFile("/home/iwsatlas1/emberger/Geant4/Current/4piSensitiveDetector_mod/B4-build/pi0_kinematics.root");
//         //f.ls();
//         Double_t energy;
//         Double_t cos_theta;
//
//         TH2D * hist1 = (TH2D*)f->Get("pi0_energyAngle");
//
//         TH1D * sample1= new TH1D("sample1", "sample1", 100, 0, 1);
//         TH1D * sample2= new TH1D("sample2", "sample2", 200, -1, 1);
//
//         TRandom3 * rndGen = new TRandom3();
//
//         TVector3 BeamDirection(0.,0.,1.);
//         BeamDirection=BeamDirection.Unit();
//
//         TH1D * check = new TH1D("Check", "check", 200, -1, 1);
//         TH2D * fd = new TH2D("fd", "fd", 200,-2,2,200,-2,2);
//         for(Int_t i=0; i<1000000; i++) {
//                 hist1->GetRandom2(energy, cos_theta);
//                 sample1->Fill(energy);
//                 sample2->Fill(cos_theta);
//
//
//                 if(cos_theta>0) {
//                         Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t radius = TMath::Tan(theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
//                         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
//                         fd->Fill(x,y);
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), BeamDirection.z());
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//                 else if(cos_theta==0) {
//                         //Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*1;
//                         Double_t y=TMath::Sin(a)*1;
//                         fd->Fill(x,y);
//
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), 0);
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//                 else if(cos_theta<0) {
//                         Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t radius = TMath::Tan(theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
//                         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
//                         fd->Fill(x,y);
//
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), -BeamDirection.z());
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//
//         }
//
//         dummy->cd(1);
//         sample1->Draw();
//         dummy->cd(2);
//         sample2->Draw();
//         dummy->cd(3);
//         check->Draw();
//         dummy->cd(4);
//         fd->Draw();
//
//
// }



void TROOTAnalysis::PrintERes(){

        TCanvas * res = new TCanvas("EnergyResolution");
        //TCanvas * gap = new TCanvas("GapEnergy", "GapEnergy");

        //Double_t y[12]={3.945/21.74,6.73/64.4,8.89/107.2,10.56/150,12.18/193,13.38/235.8,14.5/278.4,15.96/321.4, 16.57/364.1,17.32/407.2,22.07/642.2 };
        Double_t x[2]={0.1,1.0};

        Double_t y_pbabs[2]={0.108, 0.03};
        Double_t y_pbgap[2]={0.131, 0.033};
        Double_t y_cuabs[2]={0.115,0.035};
        Double_t y_cugap[2]={0.148,0.041 };
        Double_t y_cdr[2]={0.16, 0.0516};

        TGraph * re1 = new TGraph(2, x, y_pbabs);
        re1->GetYaxis()->SetRangeUser(0.02,0.18);
        re1->GetXaxis()->SetTitle("Energy [GeV]");
        re1->GetYaxis()->SetTitle("Energy Resolution");
        re1->GetYaxis()->SetTitleOffset(1.3);
        re1->SetTitle("EnergyResolution");

        TGraph * re2= new TGraph(2,x,y_pbgap);
        TGraph * re3= new TGraph(2,x,y_cuabs);
        TGraph * re4= new TGraph(2,x,y_cugap);
        TGraph * re5= new TGraph(2,x,y_cdr);

        re1->SetMarkerColor(kGray+2);
        re1->SetMarkerStyle(20);
        re1->SetMarkerSize(2);

        re2->SetMarkerColor(kGray+2);
        re2->SetMarkerStyle(24);
        re2->SetMarkerSize(2);

        re3->SetMarkerColor(kOrange+7);
        re3->SetMarkerStyle(20);
        re3->SetMarkerSize(2);

        re4->SetMarkerColor(kOrange+7);
        re4->SetMarkerStyle(24);
        re4->SetMarkerSize(2);

        re5->SetMarkerColor(kGray+2);
        re5->SetMarkerStyle(34);
        re5->SetMarkerSize(2);

        re1->Draw("AP");
        re2->Draw("P");
        re3->Draw("P");
        re4->Draw("P");
        re5->Draw("P");

        auto * legend2 = new TLegend(0.1,0.7,0.48,0.9);
        legend2->AddEntry(re1, "1mmPb_AbsorberFirst", "p");
        legend2->AddEntry(re2, "1mmPb_ScintillatorFirst", "p");
        legend2->AddEntry(re5, "1.75mmPb_AbsorberFirst_CDR", "p");
        legend2->AddEntry(re3, "2mmCu_AbsorberFirst", "p");
        legend2->AddEntry(re4, "2mmCu_ScintillatorFirst", "p");

        legend2->Draw();

        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.root");
        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.C");
        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.pdf");


        // TImage * img1 = TImage::Create();
        // img1->FromPad(res);
        // //img1->Scale(1000, 1000);
        // img1->WriteImage("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/GammaEnergyandSlopeScan_Analysis/ERes.png");
        // delete img1;
        //
        // Double_t Egun[12]={50,100,200,300,400,500,700,1000,2500,5000,7500,10000};
        // Double_t Egap[12]={19.25,38.47,77.51,115.9,154.2,192.9,268.8,384,956.2,1905,2836,3779};
        //
        // Double_t errx[12]={0,0,0,0,0,0,0,0,0,0,0,0};
        // Double_t erry[12]={4.55, 5.942, 9.056, 12.06, 11.59, 13.71, 17.98, 20.25, 53.25, 69.23, 94.6, 146.7};
        //
        // TF1 * Efficiency = new TF1("Efficiency", "[0]*x+[1]");
        //
        // TGraphErrors * gap1 = new TGraphErrors(12, Egun, Egap, errx, erry );
        // gap1->SetTitle("Efficiency");
        // gap1->Fit(Efficiency);
        // gap1->GetXaxis()->SetTitle("Gun Energy[MeV]");
        // gap1->GetYaxis()->SetTitle("Gap Energy[MeV]");
        // gap1->GetYaxis()->SetTitleOffset(1.3);
        // gStyle->SetOptFit();
        //
        // gap->cd(0);
        // gap1->Draw("A*");
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

Bool_t TROOTAnalysis::PCAEvent(Int_t event){

        Int_t nHits = Cevent->NHits();

        std::vector<TVector3> eigenVecs;

        std::vector<Double_t> showerCOG(3, 0.);

        //calculate the average shower position

        // TVector3 dummyvec(0.,0.,0.);
        // std::vector<TVector3> sample;
        // nHits=1000;
        // for(Int_t i=0;i<nHits; i++) {
        //         dummyvec.SetX(dummyvec.X()+0.1);
        //         dummyvec.SetY(dummyvec.Y()+rndGen->Gaus(0.0, 0.2));
        //         dummyvec.SetZ(dummyvec.Z()+rndGen->Gaus(0.0, 0.2));
        //
        //         sample.push_back(dummyvec);
        // }

        Double_t eges=0.;
        for(Int_t i=0; i<nHits; i++) {

                auto hitnormal=TransformCoordinates(Cevent->Hit(i)->X(), Cevent->Hit(i)->Y(), Cevent->Hit(i)->Z(), Cevent->Hit(i)->CalorimeterPart());
                Double_t edep=Cevent->Hit(i)->EnergyDeposit();

                std::tuple<Double_t, Double_t, Double_t> hit=std::make_tuple(std::get<2>(hitnormal), std::get<1>(hitnormal), std::get<0>(hitnormal));


                showerCOG.at(0)+=(std::get<0>(hit)*edep);
                showerCOG.at(1)+=(std::get<1>(hit)*edep);
                showerCOG.at(2)+=(std::get<2>(hit)*edep);
                eges+=edep;

                // showerCOG.at(0)+=sample[i].X();
                // showerCOG.at(1)+=sample[i].Y();
                // showerCOG.at(2)+=sample[i].Z();

        }

        showerCOG.at(0)/=eges;
        showerCOG.at(1)/=eges;
        showerCOG.at(2)/=eges;

        // showerCOG.at(0)/=nHits;
        // showerCOG.at(1)/=nHits;
        // showerCOG.at(2)/=nHits;
        //
        // showerCOG.at(0)/=1000;
        // showerCOG.at(1)/=1000;
        // showerCOG.at(2)/=1000;


        Double_t xi2 = 0.;
        Double_t xiyi = 0.;
        Double_t xizi = 0.;
        Double_t yi2 = 0.;
        Double_t yizi = 0.;
        Double_t zi2 = 0.;
        Double_t weightSum = 0.;

        //Get symmetric matrix


        for(Int_t i=0; i<nHits; i++) {

                auto hitnormal=TransformCoordinates(Cevent->Hit(i)->X(), Cevent->Hit(i)->Y(), Cevent->Hit(i)->Z(), Cevent->Hit(i)->CalorimeterPart());

                std::tuple<Double_t, Double_t, Double_t> hit=std::make_tuple(std::get<2>(hitnormal), std::get<1>(hitnormal), std::get<0>(hitnormal));

                Double_t weight=Cevent->Hit(i)->EnergyDeposit()*Cevent->Hit(i)->RadiusWeight();
                //Double_t weight= 1.;

                Double_t x=(std::get<0>(hit)-showerCOG.at(0))*weight;
                Double_t y=(std::get<1>(hit)-showerCOG.at(1))*weight;
                Double_t z=(std::get<2>(hit)-showerCOG.at(2))*weight;

                // Double_t x=(sample[i].X()-showerCOG.at(0))*weight;
                // Double_t y=(sample[i].Y()-showerCOG.at(1))*weight;
                // Double_t z=(sample[i].Z()-showerCOG.at(2))*weight;

                weightSum+=weight * weight;
                xi2 += x * x;
                xiyi += x * y;
                xizi += x * z;
                yi2 += y * y;
                yizi += y * z;
                zi2 += z * z;
        }

        Eigen::Matrix3f mat;

        mat <<  xi2,    xiyi,   xizi,
                xiyi,   yi2,    yizi,
                xizi,   yizi,   zi2;

        mat /= weightSum;


        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> decomp(mat);

        if (decomp.info() == Eigen::ComputationInfo::Success) {
                using eigenValColPair = std::pair<float, size_t>;
                std::vector<eigenValColPair> eigenValColVec = {
                        eigenValColPair(decomp.eigenvalues()(0), 0),
                        eigenValColPair(decomp.eigenvalues()(1), 1),
                        eigenValColPair(decomp.eigenvalues()(2), 2)
                };

                std::sort(eigenValColVec.begin(), eigenValColVec.end(),[](const eigenValColPair& left, const eigenValColPair& right){
                        return left.first > right.first;
                });

                Eigen::Matrix3f eigenVectors(decomp.eigenvectors());
                for (const auto &pair : eigenValColVec) {
                        TVector3 temp(
                                eigenVectors(0, pair.second),
                                eigenVectors(1, pair.second),
                                eigenVectors(2, pair.second));

                        eigenVecs.push_back(temp);
                }
                EstimatePhoton1.push_back(std::make_pair(TVector3(showerCOG.at(2), showerCOG.at(1), showerCOG.at(0)),
                                                         TVector3(eigenVecs[0].Z(), eigenVecs[0].Y(), eigenVecs[0].X())));
                return true;
        }
        else {
                std::cerr << "WARNING: PCA decompose failure in event "<<event <<"."<<"Number of hits: "<<nHits<<std::endl;
                return false;
        }



}


//---------------------------------------------------------------------------------------------------------------
void TROOTAnalysis::RejectOutliers(Int_t event,Double_t radius,Double_t fraction){

        // if(set==false){
        CurrentRadius=radius;
        //   set=true;
        // }
        // else if(set==true){
        //   CurrentRadius -= 1;
        // }
        std::cout<<"currentradius: "<<CurrentRadius<<std::endl;
        EcalTree->GetEntry(event);         //grab event from tree

        Int_t nhits=Cevent->NHits();

        //  if((NofRejectedHits/nhits) < fraction) {

        for(Int_t i =0; i<nhits; i++) {

                auto hitnormal=TransformCoordinates(Cevent->Hit(i)->X(), Cevent->Hit(i)->Y(), Cevent->Hit(i)->Z(), Cevent->Hit(i)->CalorimeterPart());

                TVector3 hitpos(std::get<0>(hitnormal), std::get<1>(hitnormal), std::get<2>(hitnormal));

                TVector3 pcadirection=EstimatePhoton1[event].second.Unit();

                pcadirection.Print();

                TVector3 distance_vector=(EstimatePhoton1[event].first - hitpos)-((EstimatePhoton1[event].first - hitpos)*pcadirection)*pcadirection;

                std::cout<<distance_vector.Dot(pcadirection)<<std::endl;

                Double_t distance=distance_vector.Mag();

                std::cout<<distance<<std::endl;

                Cevent->Hit(i)->SetRadiusWeight(distance);

                if(distance>CurrentRadius) {
                        std::cout<<"rejected"<<std::endl;
                        Cevent->Hit(i)->SetEnergyDeposit(0.0);
                        NofRejectedHits++;
                }
        }
        //         return false;
        // }
        // else{
        //         return true;
        // }
}

//--------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::CalcCOGPion(Int_t event){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

        // Variables for fit
        Double_t xerr=0;
        Double_t yerr=0;
        Double_t cgx=0;
        Double_t cgy=0;
        Double_t cgz=0;
        Double_t Eweight=0;

        // Int_t segment;
        //
        Double_t EGesPH1=0;
        // Double_t EGesPH2=0;
        EcalTree->GetEntry(event);         //grab event from tree
        Eges = Cevent->GapEnergy();
        Int_t cnh = Cevent->NHits();

        Double_t integral;
        Int_t nofLayerswithEdep=0;

        for(Int_t i=0; i<nofInnerLayers; i++) {         //loop over all inner layers in event
                std::string cp="InnerGapLV";

                for(Int_t j =0; j<cnh; j++) {         //loop over all hits in layer i

                        //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                        Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                        if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                //std::cout<<"fill: "<<Cevent->Hit(j)->PhotNr()<<"from photon mode:"<<phnr<<std::endl;
                                h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), edep);

                        }         //fill histogram with hits of layer i
                }
                //std::cout<<photonNR<<" : "<<std::endl;
                integral=h1->Integral();

                if(integral!=0) {         //check if layer containes energy
                        nofLayerswithEdep++;

                        Eweight=integral/Cevent->GapEnergy();  //calculate weight

                        EGesPH1+=integral;

                        cgx=h1->GetMean(1);
                        xerr=h1->GetMeanError(1);
                        //extract center of gravity and error
                        cgy=h1->GetMean(2);
                        yerr=h1->GetMeanError(2);

                        if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                        if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                        cgz=i;

                        //std::cout<<"segment: "<<std::endl;

                        std::tuple<Double_t, Double_t, Double_t> COGtransformed;

                        //h3->Fill(cgx,cgy,cgz);

                        COGtransformed=TransformCoordinates(cgx, cgy, cgz, cp);

                        auto cg=std::make_tuple(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed), xerr, yerr, Eweight);

                        //  cogs->Fill(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed));

                        //  std::cout<<"------------------------"<<std::endl;

                        coglist.push_back(cg);


                }         // end of clustering


                //ClusteredHits.clear();

                h1->Reset();

                h2->Reset();

        }
        for(Int_t i=0; i<nofOuterLayers; i++) {         //loop over all outer layers in event
                std::string cp="OuterGapLV";

                for(Int_t j =0; j<cnh; j++) {         //loop over all hits in laver i

                        //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                        Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                        if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->CalorimeterPart()==cp) {

                                //std::cout<<"fill: "<<Cevent->Hit(j)->PhotNr()<<"from photon mode:"<<phnr<<std::endl;
                                h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), edep);

                        }         //fill histogram with hits of layer i
                }
                //std::cout<<photonNR<<" : "<<std::endl;
                integral=h1->Integral();

                if(integral!=0) {         //check if layer containes energy
                        nofLayerswithEdep++;

                        Eweight=integral/Cevent->GapEnergy();  //calculate weight

                        EGesPH1+=integral;

                        cgx=h1->GetMean(1);
                        xerr=h1->GetMeanError(1);
                        //extract center of gravity and error
                        cgy=h1->GetMean(2);
                        yerr=h1->GetMeanError(2);

                        if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                        if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                        cgz=i;

                        //std::cout<<"segment: "<<std::endl;

                        std::tuple<Double_t, Double_t, Double_t> COGtransformed;

                        //h3->Fill(cgx,cgy,cgz);

                        COGtransformed=TransformCoordinates(cgx, cgy, cgz, cp);

                        auto cg=std::make_tuple(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed), xerr, yerr, Eweight);

                        //  std::cout<<"------------------------"<<std::endl;
                        //cogs->Fill(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed));

                        coglist.push_back(cg);


                }         // end of clustering



                h1->Reset();

                h2->Reset();

        }                 //end of event

        //showerCOGPhoton1.push_back(TransformCoordinates(h3->GetMean(1), h3->GetMean(2), h3->GetMean(3) /*, segment*/));
        COGCollectionPH1.push_back(coglist);
        EnergyPhoton1.push_back(std::make_pair(EGesPH1/0.52, 1 /*segment*/));
        coglist.clear();

}

//-------------------------------------------------------------------------------------------------------------------------------


std::tuple<Double_t, Double_t, Double_t > TROOTAnalysis::TransformCoordinates(Double_t x, Double_t y, Double_t z, std::string cp){

        std::tuple<Double_t,Double_t, Double_t> transformedCOGs;

        //transform Inner Calorimeter
        TVector3 coordinates;
        if(cp=="InnerGapLV") {

                Double_t offsetX=-1*calsizeXY/2+InnertiledimX/2;
                Double_t offsetY=-1*calsizeXY/2+InnertiledimY/2;
                Double_t offsetZ;

                if(Innerabsofirst) {
                        offsetZ=InnerAbsoThickness+ InnerGapThickness/2;
                        //std::cout<<"true"<<std::endl;
                }
                else{
                        offsetZ=InnerGapThickness/2;
                        //  std::cout<<"false"<<std::endl;

                }

                //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

                coordinates.SetXYZ((offsetX + x * InnertiledimX),              // X
                                   (offsetY + y * InnertiledimY),                // Y
                                   (offsetZ + z * (InnerAbsoThickness+InnerGapThickness)));     // Z
        }
        else if(cp=="OuterGapLV") {


                Double_t offsetX=-1*calsizeXY/2+OutertiledimX/2;
                Double_t offsetY=-1*calsizeXY/2+OutertiledimY/2;
                Double_t offsetZ;

                if(Outerabsofirst) {
                        offsetZ=(InnerGapThickness+InnerAbsoThickness)*nofInnerLayers+ PVesselThickness + OuterAbsoThickness + OuterGapThickness/2;
                        //std::cout<<"true"<<std::endl;
                }
                else{
                        offsetZ=(InnerGapThickness+InnerAbsoThickness)*nofInnerLayers+PVesselThickness+OuterGapThickness/2;
                        //  std::cout<<"false"<<std::endl;

                }

                //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

                coordinates.SetXYZ((offsetX + x * OutertiledimX),                        // X
                                   (offsetY + y * OutertiledimY),                          // Y
                                   (offsetZ + z * (OuterAbsoThickness+OuterGapThickness)));             // Z

        }
        //std::cout<<"Hit is in "<<cp<<std::endl;
        //coordinates.Print();

        //rotate into correct segment

        transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());


        return transformedCOGs;
}


//-------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::FitCOGsPion(Int_t event){
        gErrorIgnoreLevel = 1001;
        EcalTree->GetEntry(event);

        TVector3 ph1_orig=Cevent->MomentumPh1().Unit();

        Fcn myfcn;
        myfcn.SetMode("photon3D");

        MnUserParameters upar;
        double error_minimizer_parameters = 1e-4;


        myfcn.SetCOGs(COGCollectionPH1[event]);

        upar.Add("a_1", EstimatePhoton1[event].first.X(), error_minimizer_parameters);
        upar.Add("a_2", EstimatePhoton1[event].first.Y(), error_minimizer_parameters);
        upar.Add("a_3", EstimatePhoton1[event].first.Z(), error_minimizer_parameters);

        // upar.Fix("a_1");
        // upar.Fix("a_2");
        // upar.Fix("a_3");

        //std::cout<<"set COGPH1"<<std::endl;
        upar.Add("v_1", EstimatePhoton1[event].second.X(), error_minimizer_parameters);
        upar.Add("v_2", EstimatePhoton1[event].second.Y(), error_minimizer_parameters);
        upar.Add("v_3", EstimatePhoton1[event].second.Z(), error_minimizer_parameters);

        //  upar.SetLimits("v_3", 0.6, 1.0);
        //cout << upar << endl;

        MnMigrad migrad(myfcn, upar, 2);

        myfcn.SetCurrentEvent(event);

        //std::cout<<"event: "<<events<<std::endl;
        FunctionMinimum min = migrad();

        //std::cout<<min<<std::endl;

        MnUserParameterState uParState = min.UserState();

        //std::cout<<"attempt adding params for photon 1"<<std::endl;

        TVector3 A_1(uParState.Value("a_1"), uParState.Value("a_2"), uParState.Value("a_3"));
        TVector3 V_1(uParState.Value("v_1"), uParState.Value("v_2"), uParState.Value("v_3"));
        V_1=V_1.Unit();

        auto tp1 = std::make_pair(A_1, V_1);
        //std::cout<<"errga"<<std::endl;
        DirectionPhoton1.push_back(tp1);
        //std::cout<<"adding params for photon 1"<<std::endl;

        //TVector3 true_dir=TVector3(0.4,0.4,1.0).Unit();

        dx1->Fill(DirectionPhoton1[event].second.X());
        dy1->Fill(DirectionPhoton1[event].second.Y());
        dz1->Fill(TMath::Sqrt(DirectionPhoton1[event].second.Z()*DirectionPhoton1[event].second.Z()));

}

//---------------------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::PlotProjection(Double_t distance, Int_t event){

        Double_t k=(-distance-EstimatePhoton1[event].first.Z())/EstimatePhoton1[event].second.Z();

        TVector3 Intersect_pca(EstimatePhoton1[event].first.X()+k*EstimatePhoton1[event].second.X(),
                               EstimatePhoton1[event].first.Y()+k*EstimatePhoton1[event].second.Y(),
                               EstimatePhoton1[event].first.Z()+k*EstimatePhoton1[event].second.Z());

        projection_pca->Fill(Intersect_pca.X(), Intersect_pca.Y());

        angledifference_pca->Fill(Intersect_pca.Angle(projection_true));



        k=(-distance-DirectionPhoton1[event].first.Z())/DirectionPhoton1[event].second.Z();

        TVector3 Intersect_min(DirectionPhoton1[event].first.X()+k*DirectionPhoton1[event].second.X(),
                               DirectionPhoton1[event].first.Y()+k*DirectionPhoton1[event].second.Y(),
                               DirectionPhoton1[event].first.Z()+k*DirectionPhoton1[event].second.Z());

        projection_minimization->Fill(Intersect_min.X(), Intersect_min.Y());

        angledifference_minimization->Fill(Intersect_min.Angle(projection_true));

        // projection_correlationX->Fill(Intersect_pca.X(), Intersect_min.X());
        // projection_correlationY->Fill(Intersect_pca.Y(), Intersect_min.Y());
        //
        // projection_correlationDeltaR->Fill((Intersect_pca-projection_true).Mag(), (Intersect_min-projection_true).Mag());
        // projection_correlationDeltaRdiff->Fill((Intersect_pca-projection_true).Mag() - (Intersect_min-projection_true).Mag());
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::DrawHists(){

        dx1->GetXaxis()->SetTitle("X component of Direction Vector");
        dx1->GetYaxis()->SetTitle("Entries");

        dx2->GetXaxis()->SetTitle("X component of Direction Vector");
        dx2->GetYaxis()->SetTitle("Entries");

        dy1->GetXaxis()->SetTitle("Y component of Direction Vector");
        dy1->GetYaxis()->SetTitle("Entries");

        dy2->GetXaxis()->SetTitle("Y component of Direction Vector");
        dy2->GetYaxis()->SetTitle("Entries");

        dz1->GetXaxis()->SetTitle("Z component of Direction Vector");
        dz1->GetYaxis()->SetTitle("Entries");

        dz2->GetXaxis()->SetTitle("Z component of Direction Vector");
        dz2->GetYaxis()->SetTitle("Entries");

        D->Divide(3,2,0.01,0.01);
        gStyle->SetOptStat(2222);
        D->cd(1);
        dx1->Draw();
        D->cd(2);
        dy1->Draw();
        D->cd(3);
        dz1->Draw();
        D->cd(4);
        dx2->Draw();
        D->cd(5);
        dy2->Draw();
        D->cd(6);
        dz2->Draw();





        //projectionC->Divide(2,1,0.01,0.01);

        // projectionC->cd(1);
        // projection_pca->Draw("colz");
        // line_x->SetLineColor(kRed);
        // line_x->SetLineWidth(2);
        // line_y->SetLineColor(kRed);
        // line_y->SetLineWidth(2);
        // line_x->Draw();
        // line_y->Draw();
        //
        // projectionC->cd(2);
        // projection_minimization->Draw("colz");
        // line_x->SetLineColor(kRed);
        // line_x->SetLineWidth(2);
        // line_y->SetLineColor(kRed);
        // line_y->SetLineWidth(2);
        // line_x->Draw();
        // line_y->Draw();

        // correlationC->Divide(2,2, 0.01,0.01);
        //
        // correlationC->cd(1);
        // projection_correlationX->GetXaxis()->SetTitle("X coordinate of penetration point(orth. distance regeression)[mm]");
        // projection_correlationX->GetYaxis()->SetTitle("X coordinate of penetration point(PCA)[mm]");
        // projection_correlationX->GetYaxis()->SetTitleOffset(1.4);
        //
        // std::unique_ptr<TF1> f1(new TF1("f1", "[0]*x+[1]"));
        // //  projection_correlationX->Fit("f1");
        //
        // projection_correlationX->Draw("colz");
        //
        // correlationC->cd(2);
        // projection_correlationY->GetXaxis()->SetTitle("Y coordinate of penetration point(orth. distance regeression)[mm]");
        // projection_correlationY->GetYaxis()->SetTitle("Y coordinate of penetration point(PCA)[mm]");
        // projection_correlationY->GetYaxis()->SetTitleOffset(1.4);
        //
        // //projection_correlationY->Fit("f1");
        //
        // projection_correlationY->Draw("colz");
        //
        //
        //
        // projection_correlationDeltaR->GetXaxis()->SetTitle("Distance to true projection point(orth. distance regeression)[mm]");
        // projection_correlationDeltaR->GetYaxis()->SetTitle("Distance to true projection point(PCA)[mm]");
        // projection_correlationDeltaR->GetYaxis()->SetTitleOffset(1.4);
        //
        // correlationC->cd(3);
        // projection_correlationDeltaR->Draw("colz");
        //
        // projection_correlationDeltaRdiff->GetXaxis()->SetTitle("distance_to_true(PCA)-distance_to_true(ODR)[mm]");
        // projection_correlationDeltaRdiff->GetYaxis()->SetTitle("Entries");
        // projection_correlationDeltaRdiff->GetYaxis()->SetTitleOffset(1.4);
        //
        // correlationC->cd(4);
        // projection_correlationDeltaRdiff->Draw("colz");


        std::string e=std::to_string(true_energy);

        std::string filename="Projection_E:"+e.substr(0,e.find_last_of("."))+"MeV_";


        projectionCPCA->cd(0);
        projection_pca->Draw("colz");

        std::string save;

        save=savepath+"/ROOT/PCA"+filename+".root";
        projectionCPCA->Print(save.c_str());

        save=savepath+"/C/PCA"+filename+".C";
        projectionCPCA->Print(save.c_str());

        save=savepath+"/PDF/PCA"+filename+".pdf";
        projectionCPCA->Print(save.c_str());


        projectionCODR->cd(0);
        projection_minimization->Draw("colz");

        save=savepath+"/ROOT/ODR"+filename+".root";
        projectionCODR->Print(save.c_str());

        save=savepath+"/C/ODR"+filename+".C";
        projectionCODR->Print(save.c_str());

        save=savepath+"/PDF/ODR"+filename+".pdf";
        projectionCODR->Print(save.c_str());




        filename="Angledifference_E:"+e.substr(0,e.find_last_of("."))+"MeV_";

        angledifferenceCPCA->cd(0);
        angledifference_pca->Draw("colz");

        save=savepath+"/ROOT/PCA"+filename+".root";
        angledifferenceCPCA->Print(save.c_str());

        save=savepath+"/C/PCA"+filename+".C";
        angledifferenceCPCA->Print(save.c_str());

        save=savepath+"/PDF/PCA"+filename+".pdf";
        angledifferenceCPCA->Print(save.c_str());

        angledifferenceCODR->cd(0);
        angledifference_minimization->Draw("colz");

        save=savepath+"/ROOT/ODR"+filename+".root";
        angledifferenceCODR->Print(save.c_str());

        save=savepath+"/C/ODR"+filename+".C";
        angledifferenceCODR->Print(save.c_str());

        save=savepath+"/PDF/ODR"+filename+".pdf";
        angledifferenceCODR->Print(save.c_str());

        // TCanvas * cog= new TCanvas("cg", "cg");
        //
        // cog->cd(0);
        // //cogs->SetMarkerSize(2);
        // cogs->Draw();

        // showerdepth->cd(0);
        // showerdepthH->Draw();
//
//         Delta1->Divide(2,2,0.01,0.01);
//
//         gStyle->SetOptStat(2222);
//
//         delx->GetXaxis()->SetTitle("DeltaX[mm]");
//         dely->GetXaxis()->SetTitle("DeltaY[mm]");
//         delz->GetXaxis()->SetTitle("DeltaZ[mm]");
//
//         delx->GetYaxis()->SetTitle("Entries");
//         dely->GetYaxis()->SetTitle("Entries");
//         delz->GetYaxis()->SetTitle("Entries");
//         appdist1->GetXaxis()->SetTitle("Distance of closest Approach [mm]");
//         appdist1->GetYaxis()->SetTitle("Entries");
//
//         Delta1->cd(1);
//         delx->Draw();
//         Delta1->cd(2);
//         dely->Draw();
//         Delta1->cd(3);
//         delz->Draw();
//         Delta1->cd(4);
//         appdist1->Draw();
//
//         // Delta2->cd(0);
//         // delx->Draw();
//         //
//         // Delta3->cd(0);
//         // appdist1->Draw();
//         //
//         // Delta4->cd(0);
//         // delz->Draw();
//
//

}
