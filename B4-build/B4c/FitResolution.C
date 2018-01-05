
{
std::string folder="~/Geant4/Data/ForwardEcalWithAirGap/AngularResolution_OutlierRejection/AngReso_IT20_OT20_Ogapfirst_30Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/";
std::string file="EnergyResolution.root";

std::string path=folder+file;


TFile * f1= new TFile(path.c_str());
TCanvas * c1=(TCanvas*)f1->Get("c2p");
TGraph * g1=(TGraph*)c1->GetListOfPrimitives()->At(1);

TF1 * resolution = new TF1("resolution", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

gStyle->SetOptFit(1);
g1->Fit(resolution);

TCanvas * c2=new TCanvas();
g1->Draw("ap");


std::string filename="EnergyResolution.C";

std::string savepath=folder+filename;

c2->Print(savepath.c_str());

filename="EnergyResolution.pdf";

savepath=folder+filename;

c2->Print(savepath.c_str());
}