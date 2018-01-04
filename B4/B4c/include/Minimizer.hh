#ifndef MINIMIZER_HH
#define MINIMIZER_HH

#include <iostream>
#include <vector>
// root include
#include <TFile.h>
#include <TF1.h>
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
//minuit includes
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"
#include <Math/Vector3D.h>

using namespace ROOT::Math;
using namespace ROOT::Minuit2;
using namespace std;



class Fcn : public ROOT::Minuit2::FCNBase
{
public:
//
// /** Funktion die minimiert werden soll */
virtual double operator()(const std::vector<double> &par) const;

//
// /** Error definition: Der Fehler e auf einen Parameter x ist definiert
// * so dass f(x+-e) = f(x) + Up().
// * Für einen Chi^2 fit ist 1 das richtige um 1 sigma Fehler zu bekommen. */
virtual double Up() const;

inline void SetMode(std::string m){
        mode=m;
};
inline void SetCurrentEvent(Int_t ce){
        currentEvent=ce;
};
inline void SetSegment(Int_t s){
        seg=s;
}
inline void SetCOGs(std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > c){
        COGs=c;
        // for(int i = 0; i<COGs.size(); i++) {
        //         std::cout<<"X: "<<std::get<0>(COGs[i])<<" Y: "<<std::get<1>(COGs[i])<<" Z: "<<std::get<2>(COGs[i])<<std::endl;
        // }
};
inline void SetParamsPions(std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> >  PPions){
        ParamsPions=PPions;
};

inline void SetTrack(std::pair<TVector3, TVector3> t1, std::pair<TVector3, TVector3> t2){
        trackph1=t1;
        trackph2=t2;
}

inline void SetshowerCOG(std::tuple<Double_t,Double_t,Double_t>  c1, std::tuple<Double_t,Double_t,Double_t>  c2){
        cogph1=c1;
        cogph2=c2;
}

inline void SetDeltas(Double_t d1, Double_t d2){
        delta1=d1;
        delta2=d2;
}


inline void SetPhotonPoints( Double_t px1, Double_t py1, Double_t pz1,Double_t e1,
                             Double_t px2, Double_t py2, Double_t pz2, Double_t e2){


        l_ph1[0]=px1;
        l_ph1[1]=py1;
        l_ph1[2]=pz1;
        l_ph1[3]=e1;

        l_ph2[0]=px2;
        l_ph2[1]=py2;
        l_ph2[2]=pz2;;
        l_ph2[3]=e2;



};

inline void SetClosestApproach(Double_t cax, Double_t cay, Double_t caz){
        closestApp[0]=cax;
        closestApp[1]=cay;
        closestApp[2]=caz;

}


private:
std::tuple<Double_t,Double_t,Double_t> cogph1;
std::tuple<Double_t,Double_t,Double_t> cogph2;

std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > COGs;
std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > ParamsPions;
std::pair<TVector3, TVector3> trackph1;
std::pair<TVector3, TVector3> trackph2;


Double_t l_ph1[4];
Double_t l_ph2[4];

Double_t closestApp[4];

Double_t delta1;
Double_t delta2;



Int_t currentEvent;
Int_t seg;
std::string mode;
};



#endif  // MINIMIZER_H
