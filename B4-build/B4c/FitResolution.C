
{
std::string folder="~/Geant4/Data/ForwardEcalWithAirGap/AngularResolution_OutlierRejection/AngReso_IT20_OT20_Ogapfirst_30Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/";
std::string file="RejectionMethodAnalysis/ROOT/Res_ODR_68Quantil.root";

std::string path=folder+file;


TFile * f1= new TFile(path.c_str());
TCanvas * c1=(TCanvas*)f1->Get("c2p");
TGraph * g1=(TGraph*)c1->GetListOfPrimitives()->Last();

TF1 * resolution = new TF1("resolution", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

gStyle->SetOptFit(1);
g1->Fit(resolution);

TCanvas * c2=new TCanvas();
g1->Draw("ap");


std::string filename="AngularResolution_Titan.C";

std::string savepath=folder+filename;

c2->Print(savepath.c_str());

filename="AngularResolution_Titan.pdf";

savepath=folder+filename;

c2->Print(savepath.c_str());

std::string folder="~/Geant4/Data/ForwardEcalWithAirGap/AngularResolution_OutlierRejection/AngReso_IT20_OT20_Ogapfirst_30Inner_50Outer_lead1mm_Polystyrene5mm_20mmSteelVessel/";
std::string file="RejectionMethodAnalysis/ROOT/Res_ODR_68Quantil.root";

std::string path=folder+file;


TFile * f1= new TFile(path.c_str());
TCanvas * c1=(TCanvas*)f1->Get("c2p");
TGraph * g1=(TGraph*)c1->GetListOfPrimitives()->Last();

TF1 * resolution = new TF1("resolution", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

gStyle->SetOptFit(1);
g1->Fit(resolution);

TCanvas * c2=new TCanvas();
g1->Draw("ap");


std::string filename="AngularResolution_Steel.C";

std::string savepath=folder+filename;

c2->Print(savepath.c_str());

filename="AngularResolution_Steel.pdf";

savepath=folder+filename;

c2->Print(savepath.c_str());
}
