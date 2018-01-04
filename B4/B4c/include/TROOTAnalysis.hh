
#ifndef TROOTAnalysis_hh
#define TROOTAnalysis_hh

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"

#include "/remote/pcilc7/software/Eigen/eigen-eigen-5a0156e40feb/Eigen/Dense"

#pragma GCC diagnostic pop

#include "B4ROOTEvent.hh"

#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TBrowser.h"
#include "TColor.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TAttMarker.h"
#include "TFile.h"
#include "TList.h"
#include "TF1.h"
#include "TMath.h"
#include "TImage.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TLegend.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minimizer.hh"
#include "Minuit2/Minuit2Minimizer.h"

#include <chrono>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <memory>
#include <utility>

using namespace ROOT::Math;
using namespace ROOT::Minuit2;

class TROOTAnalysis {

public:
TROOTAnalysis();
TROOTAnalysis(std::unique_ptr<TChain> &ch, Double_t prodist);

~TROOTAnalysis();

Int_t GetNofEntries();

Bool_t PCAEvent(Int_t event);

void ApplySmearingAndCut(Double_t c, Double_t mean, Double_t sig);
void RejectOutliers(Int_t event, Double_t radius, Double_t fraction);
void PrintEdep();
void ResSummary();
void CalcAngularResolution(std::string pth, std::string title);
void PrintERes();
void PlotRMSx();
void plotEvent(Int_t pev);
void PlotProjection(Double_t distance, Int_t event);    //distance from frontface in mm
void CalcCOGPion(Int_t event);
void FitCOGsPion(Int_t event);
void DrawHists();
void PrintFitHists(Int_t minevent, Int_t maxevent);

inline void SetPath(std::string path){
        savepath=path;
        pathset=true;
}

std::tuple<Double_t, Double_t, Double_t> TransformCoordinates(Double_t x, Double_t y, Double_t z, std::string cp);

std::vector<std::pair<Double_t, Int_t> > EnergyPhoton1;
std::vector<std::pair<TVector3, TVector3> > EstimatePhoton1;
std::vector<std::pair<TVector3, TVector3> > DirectionPhoton1;
Bool_t flg;
std::unique_ptr<TTree> EcalTree;

private:

std::unique_ptr<TRandom3> rndGen;

B4ROOTEvent * Cevent;

Int_t NofRejectedHits;
Double_t CurrentRadius;
Bool_t set=false;


Int_t nofEntries;     // number of events in Tree

Double_t Eges;  //Energy in event

Double_t InnertiledimX;
Double_t InnertiledimY;
Double_t InnerAbsoThickness;
Double_t InnerGapThickness;
Double_t nofInnerLayers;

Double_t OutertiledimX;
Double_t OutertiledimY;
Double_t OuterAbsoThickness;
Double_t OuterGapThickness;
Double_t nofOuterLayers;

Double_t calsizeXY;
Double_t PVesselThickness;

Double_t histsizeX;
Double_t histsizeY;
Double_t histsizeZ;

std::string savepath;
std::string filename;

Bool_t Innerabsofirst;
Bool_t Outerabsofirst;

Bool_t pathset=false;

TVector3 gunposition;
TVector3 true_direction;
Double_t true_energy;

TVector3 projection_true;

std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > coglist;

std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > COGCollectionPH1;

std::vector<std::tuple<Double_t, Double_t, Double_t> > showerCOGPhoton1;

std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParamsGamma;

std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParams;//to be deleted, needed for old analysis

std::vector<std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > > FitParamsPions;

TH3D * hA;
TH3D * hB;

//std::unique_ptr<TH3D> h(new TH3D("ECalEvent","ECalEvent",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ));
std::unique_ptr<TH2D> h1;
std::unique_ptr<TH2D> h2;
std::unique_ptr<TH3D> h3;

//Direction difference
std::unique_ptr<TCanvas> D;
std::unique_ptr<TH1D> dx1;
std::unique_ptr<TH1D> dy1;
std::unique_ptr<TH1D> dz1;
std::unique_ptr<TH1D> dx2;
std::unique_ptr<TH1D> dy2;
std::unique_ptr<TH1D> dz2;

//hists for backprojection

std::unique_ptr<TCanvas> projectionCODR;
std::unique_ptr<TCanvas> projectionCPCA;

std::unique_ptr<TH2D> projection_pca;
std::unique_ptr<TH2D> projection_minimization;

std::unique_ptr<TCanvas> angledifferenceCPCA;
std::unique_ptr<TCanvas> angledifferenceCODR;

std::unique_ptr<TH1D> angledifference_pca;
std::unique_ptr<TH1D> angledifference_minimization;

std::unique_ptr<TCanvas> correlationC;
std::unique_ptr<TH2D> projection_correlationX;
std::unique_ptr<TH2D> projection_correlationY;


std::unique_ptr<TH2D> projection_correlationDeltaR;
std::unique_ptr<TH1D> projection_correlationDeltaRdiff;


std::unique_ptr<TLine> line_x;
std::unique_ptr<TLine> line_y;

//TH3D * cogs= new TH3D("cogs", "cogs",200,-400,0,200,-400,0,500,0,500);

// TCanvas * showerdepth = new TCanvas("showerDepthC", "showerDepthC");
// TH1D * showerdepthH = new TH1D("showerDepth", "showerDepth", 100,0,100);


};
#endif
