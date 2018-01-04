#include "PlotQuantiles.hh"


void PlotQuantiles::PlotPCA(std::string path, std::string tag){

        TCanvas * ca1p=new TCanvas("c2p","c2p");

        Int_t Estart=50;
        Double_t Xp[11];
        Double_t RmsXp[11];

        std::string fileName;

        for(Int_t i = 0; i<10; i++) {

                fileName=path+"/PCAAngledifference_E:"+std::to_string(Estart)+"MeV_.root";

                // std::cout<<t<<std::endl;
                // std::string it=t.path().string();
                // std::string filename=it.substr(it.find_last_of("/\\")+1);
                // std::cout<<filename<<std::endl;
                TFile * f1p=new TFile(fileName.c_str());

                TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCPCA");
                TH1D * hi1p = (TH1D*)c1p->GetListOfPrimitives()->Last();

                Double_t quantilePos;
                Double_t quantile=0.68;

                hi1p->GetQuantiles(1,&quantilePos,&quantile);
                std::cout<<quantilePos<<std::endl;
                RmsXp[i]=quantilePos;

                Xp[i]=Estart;

                Estart+=100;
        }

        Estart=1500;
        fileName=path+"/PCAAngledifference_E:"+std::to_string(Estart)+"MeV_.root";

        TFile * f1p=new TFile(fileName.c_str());
        TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCPCA");
        TH1D * hi1p = (TH1D*)c1p->GetListOfPrimitives()->Last();

        Double_t quantilePos;
        Double_t quantile=0.68;
        hi1p->GetQuantiles(1,&quantilePos,&quantile);
        std::cout<<quantilePos<<std::endl;

        RmsXp[10]=quantilePos;

        Xp[10]=Estart;

        TGraph * g1p = new TGraph(11, Xp, RmsXp);

        std::string plottitle=tag+"_68Quantile_PCA";

        g1p->SetTitle(plottitle.c_str());

//TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");

//g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1p->GetXaxis()->SetTitle("Energy[MeV]");
        g1p->GetYaxis()->SetTitle("AngularResolution");

        gStyle->SetOptFit();

        ca1p->cd(0);
        g1p->Draw("A*");

        std::string pathp=path+"/"+tag+"_PCA_68Quantil.root";
        ca1p->Print(pathp.c_str());

        pathp=path+"/"+tag+"_PCA_68Quantil.C";
        ca1p->Print(pathp.c_str());

        pathp=path+"/"+tag+"_PCA_68Quantil.pdf";
        ca1p->Print(pathp.c_str());

}

void PlotQuantiles::PlotODR(std::string path, std::string tag){

        TCanvas * ca1p=new TCanvas("c2p","c2p");

        Int_t Estart=50;
        Double_t Xp[11];
        Double_t RmsXp[11];

        std::string fileName;

        for(Int_t i = 0; i<10; i++) {

                fileName=path+"/ODRAngledifference_E:"+std::to_string(Estart)+"MeV_.root";

                // std::cout<<t<<std::endl;
                // std::string it=t.path().string();
                // std::string filename=it.substr(it.find_last_of("/\\")+1);
                // std::cout<<filename<<std::endl;
                TFile * f1p=new TFile(fileName.c_str());

                TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCODR");
                TH1D * hi1p = (TH1D*)c1p->GetListOfPrimitives()->Last();

                Double_t quantilePos;
                Double_t quantile=0.68;

                hi1p->GetQuantiles(1,&quantilePos,&quantile);
                std::cout<<quantilePos<<std::endl;
                RmsXp[i]=quantilePos;

                Xp[i]=Estart;

                Estart+=100;
        }

        Estart=1500;
        fileName=path+"/ODRAngledifference_E:"+std::to_string(Estart)+"MeV_.root";

        TFile * f1p=new TFile(fileName.c_str());
        TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCODR");
        TH1D * hi1p = (TH1D*)c1p->GetListOfPrimitives()->Last();

        Double_t quantilePos;
        Double_t quantile=0.68;
        hi1p->GetQuantiles(1,&quantilePos,&quantile);
        std::cout<<quantilePos<<std::endl;

        RmsXp[10]=quantilePos;

        Xp[10]=Estart;

        TGraph * g1p = new TGraph(11, Xp, RmsXp);

        std::string plottitle=tag+"_68Quantile_ODR";

        g1p->SetTitle(plottitle.c_str());

//TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");

//g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1p->GetXaxis()->SetTitle("Energy[MeV]");
        g1p->GetYaxis()->SetTitle("AngularResolution");

        gStyle->SetOptFit();

        ca1p->cd(0);
        g1p->Draw("A*");

        std::string pathp=path+"/"+tag+"_ODR_68Quantil.root";
        ca1p->Print(pathp.c_str());

        pathp=path+"/"+tag+"_ODR_68Quantil.C";
        ca1p->Print(pathp.c_str());

        pathp=path+"/"+tag+"_ODR_68Quantil.pdf";
        ca1p->Print(pathp.c_str());

}
