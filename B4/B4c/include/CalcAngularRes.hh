
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
#include <cstddef>

// #include <boost/filesystem.hpp>
//
// using namespace boost::filesystem;
// struct recursive_directory_range
// {
//     typedef recursive_directory_iterator iterator;
//     recursive_directory_range(path p) : p_(p) {}
//
//     iterator begin() { return recursive_directory_iterator(p_); }
//     iterator end() { return recursive_directory_iterator(); }
//
//     path p_;
// };



class AngularResolution {

public:
void CalcAngularResolution(std::string pth, std::string title);
void EresAndShowerprofile(std::string pth);




private:

Int_t nofEvents;

};
