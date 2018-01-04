// Main function for analyzing the TTree

#include <iostream>
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"
#include "TROOTAnalysis.hh"
//#include "PrincipalComponentAnalysis.h"
#include <stdlib.h>
#include <memory>
#include "TMath.h"


int main(int argc, char * argv[]) {

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        std::unique_ptr<TApplication> app(new TApplication("app", 0, 0, 0));
        //auto app=make_unique<TApplication>("app", 0, 0, 0);
        //TApplication* app = new TApplication("app", 0, 0, 0);
//
//    gStyle->SetCanvasBorderMode(0);
//    gStyle->SetPadBorderMode(0);
//
//    gStyle->SetPalette(1);
//    gStyle->SetOptStat(1);
//    gStyle->SetOptFit(11);
//    //gStyle->SetOptTitle(0);
//    gStyle->SetStatBorderSize(1);
//    gStyle->SetStatColor(10);
//    gStyle->SetCanvasColor(10);
//    gStyle->SetPadLeftMargin(0.16);
//    gStyle->SetPadBottomMargin(0.16);
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
//    gStyle->SetOptTitle(0);
//    gStyle->SetTitleSize(0.048,"xy");
//    gStyle->SetLabelSize(0.04,"xy");
//    gStyle->SetTitleOffset(1.3,"x");
//    gStyle->SetTitleOffset(1.3,"y");


        Double_t cut=std::stod(argv[1]);
        Double_t distance=std::stod(argv[2]);
        std::string treepath= argv[3];


        std::unique_ptr<TChain> ch1(new TChain("eventTree"));

        //std::cout<<argv[3]<<std::endl;
        ch1->Add(treepath.c_str());
        ch1->Draw("");
        //std::cout<<"hello"<<std::endl;
        TROOTAnalysis A(ch1,distance);

        A.SetPath(argv[4]);

        Double_t entries=A.GetNofEntries();

        Int_t eventA=0;

        Int_t progress=0;

        std::unique_ptr<TCanvas> cdir(new TCanvas("PCAReconstruction"));

        cdir->Divide(3,1,0.01,0.01);

        std::unique_ptr<TH1D> hx(new TH1D("xdirection", "xdirection", 1000,-3,3));
        std::unique_ptr<TH1D> hy(new TH1D("ydirection", "ydirection", 1000,-3,3));
        std::unique_ptr<TH1D> hz(new TH1D("zdirection", "zxdirection", 1000,-3,3));

        hx->GetXaxis()->SetTitle("X component of Direction Vector");
        hy->GetXaxis()->SetTitle("Y component of Direction Vector");
        hz->GetXaxis()->SetTitle("Z component of Direction Vector");

        hx->GetYaxis()->SetTitle("Entries");
        hy->GetYaxis()->SetTitle("Entries");
        hz->GetYaxis()->SetTitle("Entries");




        //entries=100000;
        for(Int_t i=0; i<entries; i++) {

                A.EcalTree->GetEntry(i);

                A.ApplySmearingAndCut(cut, 0,0.1);
                //A.PrintEdep();
                if(A.PCAEvent(i)) {

                        // Bool_t donerejection=false;
                        //
                        // while(!donerejection) {
                                //std::cout<<"iteration"<<std::endl;
                              A.RejectOutliers(i, 100, 0.2);
                                //
                                // if(!donerejection) {

                                        A.PCAEvent(i);
                        //
                        //         }
                        // }

                        A.CalcCOGPion(eventA);
                        A.FitCOGsPion(eventA);

                        A.PlotProjection(distance, eventA);

                        eventA++;

                }

                progress = ((Double_t)i /entries) * 100;


                if ((Int_t)progress % 5 == 0)
                {
                        std::cout << "\r" << std::string(progress/5, '|') << progress << "%";
                        std::cout.flush();
                }
        }



        // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // std::cout << "Computing took "
        //           << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        //           <<" milliseconds"<<std::endl;


        // TF1 *f1 = new TF1("f1","gaus",-1,1);
        // f1->SetRange(-1,1);
        // gStyle->SetOptFit(111111111);
        //
        // cdir->cd(1);
        // //hx->Fit(f1, "L", "",-1,1);
        // hx->Draw();
        //
        // cdir->cd(2);
        // //hy->Fit(f1, "L", "", -1,1);
        // hy->Draw();
        //
        // cdir->cd(3);
        // //hz->Fit(f1, "L", "", -1,1);
        // hz->Draw();

        A.DrawHists();

        // end = std::chrono::steady_clock::now();
        // std::cout << "Computing took "
        //           << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        //           <<" milliseconds"<<std::endl;
        //app->Run();
        return 0;
}
