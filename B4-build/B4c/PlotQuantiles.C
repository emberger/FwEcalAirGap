{

TCanvas * ca1p=new TCanvas("c2p","c2p");

Int_t Estart=50;
Double_t Xp[11];
Double_t RmsXp[11];

//std::string filename=pth.substr(pth.find_last_of("/\\")+1);

for(Int_t i = 0; i<10; i++) {

        fileName="/home/iwsatlas1/emberger/Geant4/Current/ForwardEcalWithPressureVessel/B4-build/B4c/AngReso_IT20_OT20_25Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/Analysis/ROOT/PCAAngledifference_E:"+std::to_string(Estart)+"MeV_.root";
        // std::cout<<t<<std::endl;
        // std::string it=t.path().string();
        // std::string filename=it.substr(it.find_last_of("/\\")+1);
        // std::cout<<filename<<std::endl;
        TFile * f1p=new TFile(fileName.c_str());

        TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCPCA");
        TH1D * hi1p = (TH1D*)c1p->GetPrimitive("AngleDiffPCA");

        Double_t quantilePos;

        hi1p->GetQuantiles(1,0.68,quantilePos);

        RmsXp[i]=quantilePos;


        Xp[i]=Estart;

        Estart+=100;
}

Estart=1500;
fileName="/home/iwsatlas1/emberger/Geant4/Current/ForwardEcalWithPressureVessel/B4-build/B4c/AngReso_IT20_OT20_25Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/Analysis/ROOT/PCAAngledifference_E:"+std::to_string(Estart)+"MeV_.root";

TFile * f1p=new TFile(fileName.c_str());
TCanvas * c1p =(TCanvas*)f1p->Get("angledifferenceCPCA");
TH2D * hi1p = (TH2D*)c1p->GetPrimitive("AngleDiffPCA");

Double_t quantilePos;

hi1p->GetQuantiles(1,0.68,quantilePos);

RmsXp[i]=quantilePos;

Xp[10]=Estart;

TGraph * g1p = new TGraph(11, X, RmsXp);


plottitle=title+"_68Quantile_PCA_X";

g1p->SetTitle(plottitle.c_str());


//TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");




//g1->Fit(rms_x);

gStyle->SetOptFit(11111111);
g1p->GetXaxis()->SetTitle("Energy[MeV]");
g1p->GetYaxis()->SetTitle("AngularResolution");



gStyle->SetOptFit();

ca1p->cd(0);
g1p->Draw("A*");
std::string title="Res";
std::string pathp=title+"PCA_X_68Quantil.root";
ca1p->Print(pathp.c_str());

pathp=title+"PCA_X_68Quantil.C";
ca1p->Print(pathp.c_str());

pathp=title+"PCA_X_68Quantil.pdf";
ca1p->Print(pathp.c_str());

}
