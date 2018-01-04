#ifndef B4ROOTEvent_hh
#define B4ROOTEvent_hh

// Header file for Pi0 Event class and sub-classes with specific information
// create dictionary with rootcint RootClasses_dict.cpp -c include/B4ROOTEvent.hh include/LinkDef.h


#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"
#include <iostream>

//#include "globals.hh"

class B4ROOTHit : public TObject {

private:
Double_t m_X;
Double_t m_Y;
Double_t m_Z;
std::string m_calorimeterPart;

Double_t m_EnergyDeposit;

Double_t m_radiusWeight;


public:
B4ROOTHit();
B4ROOTHit(const B4ROOTHit& orig);
virtual ~B4ROOTHit();
void Clear(const Option_t* option ="") {
};
// Getters:
Double_t X() const {
								return m_X;
}
Double_t Y() const {
								return m_Y;
}
Double_t Z() const {
								return m_Z;
}

std::string CalorimeterPart() const {
								return m_calorimeterPart;
}
Double_t EnergyDeposit() const {
								return m_EnergyDeposit;
}

Double_t RadiusWeight() const {
								return m_radiusWeight;
}

// Setters:

void SetCoordinates(Double_t x, Double_t y, Double_t z) {
								m_X = x; m_Y = y; m_Z = z;
}

void SetCalorimeterPart(std::string s){
								m_calorimeterPart=s;
}
void SetEnergyDeposit(Double_t dep) {
								m_EnergyDeposit = dep;
}

void SetRadiusWeight(Double_t rad){

								m_radiusWeight=1/rad;

}

ClassDef(B4ROOTHit,1)
};




class B4ROOTEvent : public TObject {

private:
Int_t m_EventNo;

Double_t m_GapEnergy;

Int_t m_NHits;

Double_t m_InnerAbsoThickness;
Double_t m_InnerGapThickness;
Double_t m_InnerLayerno;
Double_t m_InnertilesizeX;
Double_t m_InnertilesizeY;
Bool_t Innerabsofirst;

Double_t m_OuterAbsoThickness;
Double_t m_OuterGapThickness;
Double_t m_OuterLayerno;
Double_t m_OutertilesizeX;
Double_t m_OutertilesizeY;
Bool_t Outerabsofirst;

//Double_t m_airgapThickness

Double_t m_calsizeXY;
Double_t m_PVesselThickness;

TVector3 m_MomentumPh1;
Double_t m_EnergyPrimary;

TVector3 m_GunPos;


TClonesArray *m_Hits;  //->
static TClonesArray *aHits;

public:
B4ROOTEvent();
B4ROOTEvent(const B4ROOTEvent& origEvent);
// B4ROOTEvent &operator=(B4ROOTEvent &other){
//         std::cout<<"=called"<<std::endl;
//         new(this)B4ROOTEvent(other);
//         return *this;
//
// }
virtual ~B4ROOTEvent();
void Clear(const Option_t* option ="");
static void   Reset(Option_t *option ="");



//Getters
Int_t EventNo() const {
								return m_EventNo;
}
Double_t GapEnergy() const {
								return m_GapEnergy;
}

Double_t InnerAbsoThickness(){
								return m_InnerAbsoThickness;
}
Double_t InnerGapThickness(){
								return m_InnerGapThickness;
}
Double_t NumberOfInnerLayers(){
								return m_InnerLayerno;
}
Double_t InnerTilesizeX(){
								return m_InnertilesizeX;
}
Double_t InnerTilesizeY(){
								return m_InnertilesizeY;
}
Bool_t GetInnerAbsFirst(){
								return Innerabsofirst;
}


Double_t OuterAbsoThickness(){
								return m_OuterAbsoThickness;
}
Double_t OuterGapThickness(){
								return m_OuterGapThickness;
}
Double_t NumberOfOuterLayers(){
								return m_OuterLayerno;
}
Double_t OuterTilesizeX(){
								return m_OutertilesizeX;
}
Double_t OuterTilesizeY(){
								return m_OutertilesizeY;
}
Bool_t GetOuterAbsFirst(){
								return Outerabsofirst;
}

// Double_t airgapThickness(){
// 								return m_airgapThickness;
// }

Double_t calsizeXY(){
								return m_calsizeXY;
}

Double_t PVesselThickness(){
								return m_PVesselThickness;
}

Double_t EnergyPrimary(){
								return m_EnergyPrimary;
}

TVector3 MomentumPh1(){
								return m_MomentumPh1;
}

TVector3 GunPos(){
								return m_GunPos;
}


Int_t NHits() const {
								return m_NHits;
}

TClonesArray* Hits() const {
								return m_Hits;
}

B4ROOTHit* Hit(Int_t i) {
								return (B4ROOTHit*)m_Hits->At(i);
}

//Setters
void SetEventNo(Int_t evN) {
								m_EventNo = evN;
}
void SetGapEnergy(Double_t en) {
								m_GapEnergy = en;
}


void SetInnerAbsoThickness(Double_t at){
								m_InnerAbsoThickness=at;
}
void SetInnerGapThickness(Double_t gt){
								m_InnerGapThickness=gt;
}
void SetNumberOfInnerLayers(Double_t nol){
								m_InnerLayerno=nol;
}
void SetInnerTilesizeX(Double_t tsX){
								m_InnertilesizeX=tsX;
}
void SetInnerTilesizeY(Double_t tsY){
								m_InnertilesizeY=tsY;
}

void SetInnerAbsFirst(Bool_t cfg){
								Innerabsofirst=cfg;
}


void SetOuterAbsoThickness(Double_t at){
								m_OuterAbsoThickness=at;
}
void SetOuterGapThickness(Double_t gt){
								m_OuterGapThickness=gt;
}
void SetNumberOfOuterLayers(Double_t nol){
								m_OuterLayerno=nol;
}
void SetOuterTilesizeX(Double_t tsX){
								m_OutertilesizeX=tsX;
}
void SetOuterTilesizeY(Double_t tsY){
								m_OutertilesizeY=tsY;
}
void SetOuterAbsFirst(Bool_t cfg){
								Outerabsofirst=cfg;
}

// void SetairgapThickness(Double_t agt)
// {
// 								m_airgapThickness=agt;
// }

void SetPVesselThickness(Double_t pt){
								m_PVesselThickness=pt;
}

void SetcalsizeXY(Double_t cxy){
								m_calsizeXY=cxy;
}
void SetEnergyPrimary(Double_t ep){
								m_EnergyPrimary=ep;
}
void SetMomentumPh1(Double_t mphx, Double_t mphy, Double_t mphz){
								m_MomentumPh1.SetXYZ(mphx, mphy,mphz);
}
void SetGunPos(TVector3 gp){
								m_GunPos=gp;
}

B4ROOTHit* AddHit(B4ROOTHit& cand);

ClassDef(B4ROOTEvent,1)
};
#endif
