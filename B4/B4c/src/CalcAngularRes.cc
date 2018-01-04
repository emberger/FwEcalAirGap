#include "CalcAngularRes.hh"





void AngularResolution::CalcAngularResolution(std::string pth, std::string title){

        TCanvas * ca1=new TCanvas("c1","c1");
        TCanvas * ca2=new TCanvas("c2","c2");

        std::string fileName;
        Int_t Estart=50;
        Double_t X[11];
        Double_t RmsX[11];
        Double_t RmsY[11];


        //std::string filename=pth.substr(pth.find_last_of("/\\")+1);

        for(Int_t i = 0; i<10; i++) {

                fileName=pth+"/ODRProjection_E:"+std::to_string(Estart)+"MeV_.root";
                // std::cout<<t<<std::endl;
                // std::string it=t.path().string();
                // std::string filename=it.substr(it.find_last_of("/\\")+1);
                // std::cout<<filename<<std::endl;
                TFile * f1=new TFile(fileName.c_str());

                TCanvas * c1 =(TCanvas*)f1->Get("ProjectionCODR");
                TH2D * hi1 = (TH2D*)c1->GetPrimitive("Projection_minimization");

                TH1D * h1=hi1->ProjectionX();
                h1->Rebin((Int_t)h1->GetRMS()/2);
                TH1D * h2=hi1->ProjectionY();
                h2->Rebin((Int_t)h2->GetRMS()/2);


                TF1 * func1= new TF1("gausfunc1", "gaus",h1->GetMean()-h1->GetRMS(),h1->GetMean()+h1->GetRMS());

                h1->Fit(func1, "R");
                std::cout<<func1->GetParameter(0)<<":"<<func1->GetParameter(1)<<":"<<func1->GetParameter(2)<<std::endl;
                RmsX[i]=func1->GetParameter(2)/1150;

                TF1 * func2= new TF1("gausfunc2", "gaus",h1->GetMean()-h1->GetRMS(),h1->GetMean()+h1->GetRMS());

                h2->Fit(func2, "R");
                std::cout<<func2->GetParameter(0)<<":"<<func2->GetParameter(1)<<":"<<func2->GetParameter(2)<<std::endl;
                RmsY[i]=func2->GetParameter(2)/1150;

                // TCanvas dummyC;
                //
                // Double_t rmsx=hi1->GetRMS(1);
                //
                // Double_t rmsy=hi1->GetRMS(2);

                X[i]=Estart;



                Estart+=100;
        }

        Estart=1500;
        fileName=pth+"/ODRProjection_E:"+std::to_string(Estart)+"MeV_.root";

        TFile * f1=new TFile(fileName.c_str());
        TCanvas * c1 =(TCanvas*)f1->Get("ProjectionCODR");
        TH2D * hi1 = (TH2D*)c1->GetPrimitive("Projection_minimization");

        TH1D * h1=hi1->ProjectionX();
        h1->Rebin((Int_t)h1->GetRMS()/2);
        TH1D * h2=hi1->ProjectionY();
        h2->Rebin((Int_t)h2->GetRMS()/2);

        TF1 * func1= new TF1("gausfunc1", "gaus",h1->GetMean()-h1->GetRMS(),h1->GetMean()+h1->GetRMS());
        h1->Fit(func1, "R");
        std::cout<<func1->GetParameter(0)<<":"<<func1->GetParameter(1)<<":"<<func1->GetParameter(2)<<std::endl;
        RmsX[10]=func1->GetParameter(2)/1150;

        TF1 * func2= new TF1("gausfunc1", "gaus",h2->GetMean()-h2->GetRMS(),h2->GetMean()+h2->GetRMS());
        h2->Fit(func2, "R");
        std::cout<<func2->GetParameter(0)<<":"<<func2->GetParameter(1)<<":"<<func2->GetParameter(2)<<std::endl;
        RmsY[10]=func2->GetParameter(2)/1150;

        X[10]=Estart;

        TGraph * g1 = new TGraph(11, X, RmsX);
        TGraph * g1y = new TGraph(11, X, RmsY);

        std::string plottitle=title+"_ODR_X";

        g1->SetTitle(plottitle.c_str());

        plottitle=title+"_ODR_Y";
        g1y->SetTitle(plottitle.c_str());

        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");
        //TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");

        //g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1->GetXaxis()->SetTitle("Energy[MeV]");
        g1->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
        g1y->GetXaxis()->SetTitle("Energy[MeV]");
        g1y->GetYaxis()->SetTitle("#frac{RMS_y}{1000mm}");

        gStyle->SetOptFit();

        ca1->cd(0);
        g1->Draw("A*");

        std::string path=pth+"/"+title+"ODR_X_fit_rebin.root";
        ca1->Print(path.c_str());

        path=pth+"/"+title+"ODRX_X_fit_rebin.C";
        ca1->Print(path.c_str());

        path=pth+"/"+title+"ODR_X_fit_rebin.pdf";
        ca1->Print(path.c_str());


        ca2->cd(0);
        g1y->Draw("A*");

        path=pth+"/"+title+"ODR_Y_fit_rebin.root";
        ca2->Print(path.c_str());

        path=pth+"/"+title+"ODR_Y_fit_rebin.C";
        ca2->Print(path.c_str());

        path=pth+"/"+title+"ODR_Y_fit_rebin.pdf";
        ca2->Print(path.c_str());

        TCanvas * ca1p=new TCanvas("c1p","c1p");
        TCanvas * ca2p=new TCanvas("c2p","c2p");

        Estart=50;
        Double_t Xp[11];
        Double_t RmsXp[11];
        Double_t RmsYp[11];

        //std::string filename=pth.substr(pth.find_last_of("/\\")+1);

        for(Int_t i = 0; i<10; i++) {

                fileName=pth+"/PCAProjection_E:"+std::to_string(Estart)+"MeV_.root";
                // std::cout<<t<<std::endl;
                // std::string it=t.path().string();
                // std::string filename=it.substr(it.find_last_of("/\\")+1);
                // std::cout<<filename<<std::endl;
                TFile * f1p=new TFile(fileName.c_str());

                TCanvas * c1p =(TCanvas*)f1p->Get("ProjectionCPCA");
                TH2D * hi1p = (TH2D*)c1p->GetPrimitive("Projection_pca");

                TH1D * h1p=hi1p->ProjectionX();
                h1p->Rebin((Int_t)h1p->GetRMS()/2);
                TH1D * h2p=hi1p->ProjectionY();
                h2p->Rebin((Int_t)h2p->GetRMS()/2);

                TF1 * func1p= new TF1("gausfunc1p", "gaus",h1p->GetMean()-h1p->GetRMS(),h1p->GetMean()+h1p->GetRMS());

                h1p->Fit(func1p, "R");
                std::cout<<func1p->GetParameter(0)<<":"<<func1p->GetParameter(1)<<":"<<func1p->GetParameter(2)<<std::endl;
                RmsXp[i]=func1p->GetParameter(2)/1150;

                TF1 * func2p= new TF1("gausfunc2p", "gaus",h2p->GetMean()-h2p->GetRMS(),h2p->GetMean()+h2p->GetRMS());

                h2p->Fit(func2p, "R");
                std::cout<<func2p->GetParameter(0)<<":"<<func2p->GetParameter(1)<<":"<<func2p->GetParameter(2)<<std::endl;
                RmsYp[i]=func2p->GetParameter(2)/1150;

                Xp[i]=Estart;

                Estart+=100;
        }

        Estart=1500;
        fileName=pth+"/PCAProjection_E:"+std::to_string(Estart)+"MeV_.root";

        TFile * f1p=new TFile(fileName.c_str());
        TCanvas * c1p =(TCanvas*)f1p->Get("ProjectionCPCA");
        TH2D * hi1p = (TH2D*)c1p->GetPrimitive("Projection_pca");

        TH1D * h1p=hi1p->ProjectionX();
        h1p->Rebin((Int_t)h1p->GetRMS()/2);
        TH1D * h2p=hi1p->ProjectionY();
        h2p->Rebin((Int_t)h2p->GetRMS()/2);

        TF1 * func1p= new TF1("gausfunc1p", "gaus",h1p->GetMean()-h1p->GetRMS(),h1p->GetMean()+h1p->GetRMS());

        h1p->Fit(func1p, "R");
        std::cout<<func1p->GetParameter(0)<<":"<<func1p->GetParameter(1)<<":"<<func1p->GetParameter(2)<<std::endl;
        RmsXp[10]=func1p->GetParameter(2)/1150;

        TF1 * func2p= new TF1("gausfunc2p", "gaus",h2p->GetMean()-h2p->GetRMS(),h2p->GetMean()+h2p->GetRMS());

        h2p->Fit(func2p, "R");
        std::cout<<func2p->GetParameter(0)<<":"<<func2p->GetParameter(1)<<":"<<func2p->GetParameter(2)<<std::endl;
        RmsYp[10]=func2p->GetParameter(2)/1150;



        Xp[10]=Estart;

        TGraph * g1p = new TGraph(11, X, RmsXp);
        TGraph * g1yp = new TGraph(11, X, RmsYp);

        plottitle=title+"_PCA_X";

        g1p->SetTitle(plottitle.c_str());

        plottitle=title+"_PCA_Y";
        g1yp->SetTitle(plottitle.c_str());

        TF1 * rms_xp = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");
        //TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");




        //g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1p->GetXaxis()->SetTitle("Energy[MeV]");
        g1p->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
        g1yp->GetXaxis()->SetTitle("Energy[MeV]");
        g1yp->GetYaxis()->SetTitle("#frac{RMS_y}{1000mm}");


        gStyle->SetOptFit();

        ca1p->cd(0);
        g1p->Draw("A*");

        std::string pathp=pth+"/"+title+"PCA_X_fit_rebin.root";
        ca1p->Print(pathp.c_str());

        pathp=pth+"/"+title+"PCA_X_fit_rebin.C";
        ca1p->Print(pathp.c_str());

        pathp=pth+"/"+title+"PCA_X_fit_rebin.pdf";
        ca1p->Print(pathp.c_str());

        ca2p->cd(0);
        g1yp->Draw("A*");

        pathp=pth+"/"+title+"PCA_Y_fit_rebin.root";
        ca2p->Print(pathp.c_str());

        pathp=pth+"/"+title+"PCA_Y_fit_rebin.C";
        ca2p->Print(pathp.c_str());

        pathp=pth+"/"+title+"PCA_Y_fit_rebin.pdf";
        ca2p->Print(pathp.c_str());


}
